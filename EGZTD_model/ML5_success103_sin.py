#!/usr/bin/env python
# coding: utf-8

# 缺少的库直接安装pip/conda
import math
import codecs
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
# import xgboost as xgb  # 原生API，对应xgb_model = xgb.XGBRegressor()
import seaborn as sns
# from catboost import CatBoostRegressor
from lightgbm import LGBMRegressor, log_evaluation, early_stopping
from xgboost import XGBRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.base import BaseEstimator, RegressorMixin  # 这两个类是scikit-learn中的基类，确保自定义的模型类符合scikit-learn的接口标准
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split, KFold, GridSearchCV
from sklearn.compose import TransformedTargetRegressor
from sklearn.neighbors import KNeighborsRegressor
from tensorflow.keras import layers, models, optimizers
from tensorflow.keras.layers import Dropout
from scikeras.wrappers import KerasRegressor
from tensorflow.keras.callbacks import ReduceLROnPlateau, EarlyStopping
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l1_l2, l2


# 一个自定义的模型类，名为EarlyStoppingLGBMRegressor，目的是将早停功能与LGBMRegressor结合
# 显式声明参数：所有需要通过网格搜索调优的 LightGBM 参数必须在 __init__ 中显式定义
# 参数传递链：GridSearchCV参数 → Pipeline步骤参数 → EarlyStoppingLGBMRegressor属性 → LGBMRegressor参数。
# 通过这一机制，无论是固定参数还是待调优参数，均能正确覆盖默认值，确保模型按预期训练。
class EarlyStoppingLGBMRegressor(BaseEstimator, RegressorMixin):
    def __init__(
            self,
            # 必须显式声明需要调优的参数，之后交叉验证GridSearchCV要调试哪一个就传递哪一个
            learning_rate: float = 0.1,  # GridSearchCV 的参数验证逻辑仅识别显式定义的类属性，无法传递 model__learning_rate等参数
            max_bin: int = 255,  # 这些参数后续会被param_grid1中的数覆盖
            max_depth: int = 6,
            num_leaves: int = 31,
            # 其他需要调优的参数（早停控制参数），这部分参数是设计这个函数的主要目的，尤其是早停验证集，（就是为了这点醋才包的饺子），其余的全是用来配合的
            early_stopping_rounds: int = 6,
            eval_metric: str = 'rmse',
            # 以五折交叉验证为例，这里我们将建模数据中的80%（4折）作为训练数据，这些训练数据再挑选出10%（0.1*0.8=0.08）作为和验证集（20%，不是测试集）无关的早停验证集
            val_size: float = 0.1,
            random_state: int = 42,
            # 接收其他未显式声明的参数（如bagging_fraction等）
            **kwargs  # 用于接收未在 __init__ 中显式定义的额外参数，避免代码报错
    ):
        # 显式定义需调优的参数
        self.learning_rate = learning_rate
        self.max_bin = max_bin
        self.num_leaves = num_leaves
        self.max_depth = max_depth
        # 早停控制参数
        self.early_stopping_rounds = early_stopping_rounds  # 传入early_stopping_rounds参数
        self.eval_metric = eval_metric
        self.val_size = val_size  # 如果数据足够多的话，可以设置0.2，根据实际情况来
        self.random_state = random_state
        # 其他参数存入字典
        self.lgbm_kwargs = kwargs
        self.model = None

    def fit(self, x1, y1, eval_set=None):  # 输入训练数据x1和目标标签y1
        # 动态划分早停验证子集（从训练数据中）
        x_train_sub, x_val_sub, y_train_sub, y_val_sub = train_test_split(
            x1, y1, test_size=self.val_size, random_state=42)  # 42，123都是常用的随机种子

        # 初始化LightGBM模型（构造函数参数），将LGBM回归函数传递到EarlyStoppingLGBMRegressor中
        self.model = LGBMRegressor(
            learning_rate=self.learning_rate,
            max_bin=self.max_bin,
            max_depth=self.max_depth,
            metric=self.eval_metric,  # 早停评估指标,metric而非eval_metric
            random_state=self.random_state,
            **self.lgbm_kwargs  # 传递额外参数
        )
        # 定义LightGBM回调函数
        callbacks_lgbm = [log_evaluation(period=-1),  # 每100轮输出一次日志,-1为不输出日志
                          early_stopping(self.early_stopping_rounds)]  # 如果5轮没有提高，停止训练
        # 绑定早停验证集
        self.model.fit(
            x_train_sub, y_train_sub,  # 下行显式传入早停验证子集
            eval_set=[(x_val_sub, y_val_sub)],  # 直接传入eval_set。
            eval_metric=self.eval_metric,  # fit中采用的显示传递的指标
            callbacks=callbacks_lgbm
        )
        return self

    def predict(self, x1):
        if self.model is None:
            raise RuntimeError("The model must be trained with fit () method first")
        return self.model.predict(x1)

    # 代理 feature_importances_，此外，def __getattr__(self, name):也可以代理属性，但是更难、更强大一些
    @property
    def feature_importances_(self):
        if self.model is None:
            raise RuntimeError("The model must be trained with fit () method first")  # 防御性编程
        return self.model.feature_importances_


# XGBoost的增强版，自适应随机选取早停验证集，可配合格网搜索+交叉验证，同时传递了一大堆关键性能
class XGBoostEarlyStoppingRegressor(BaseEstimator, RegressorMixin):
    def __init__(
            self,
            # 基础参数
            learning_rate: float = 0.1,
            max_depth: int = 6,
            n_estimators: int = 100,
            # 早停控制
            early_stopping_rounds: int = 6,
            eval_metric: str = 'rmse',
            # 验证设置
            val_size: float = 0.1,
            random_state: int = 42,
            # 高级参数
            tree_method: str = 'hist',  # 内存优化选项
            # verbosity: int = 0,  # 0(silent),1(warning),2(info),3(debug)
            **xgb_kwargs  # 支持其他XGBoost原生参数
    ):
        # 基础参数，只给格网搜索中需要调优的参数即可
        self.learning_rate = learning_rate
        self.max_depth = max_depth
        self.n_estimators = n_estimators

        # 早停控制
        self.early_stopping_rounds = early_stopping_rounds
        self.eval_metric = eval_metric
        # 验证设置
        self.val_size = val_size
        self.random_state = random_state
        # 高级参数
        self.tree_method = tree_method
        # self.verbosity = verbosity
        self.xgb_kwargs = xgb_kwargs  # 保存扩展参数
        self.model = None
        # self.best_iteration_ = 0  # 记录最优迭代次数

    def fit(self, x2, y2, eval_set=None):
        # val_set - 可手动指定验证集，格式为[(X_val, y_val)], 若为None则自动划分
        if eval_set is None:  # 自动划分验证集
            x_train_sub, x_val_sub, y_train_sub, y_val_sub = train_test_split(
                x2, y2, test_size=self.val_size, random_state=self.random_state)  # 42，123都是常用的随机种子
            eval_set = [(x_val_sub, y_val_sub)]
        else:
            x_train_sub, y_train_sub = x2, y2

        # 模型初始化
        self.model = XGBRegressor(  # 这一部分优先级要高于前面self那一部分
            learning_rate=self.learning_rate,
            max_depth=self.max_depth,
            n_estimators=self.n_estimators,
            eval_metric=self.eval_metric,
            early_stopping_rounds=self.early_stopping_rounds,
            tree_method=self.tree_method,
            # verbosity=self.verbosity,
            **self.xgb_kwargs
        )

        # 执行训练（抑制警告）
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.model.fit(
                x_train_sub, y_train_sub,
                eval_set=eval_set,
                verbose=False  # 关闭内置日志
            )
        # 记录最优迭代次数
        # self.best_iteration_ = self.model.best_iteration
        return self

    def predict(self, x2):  # 使用最优迭代次数进行预测
        if self.model is None:
            raise RuntimeError("The model must be trained with fit () method first")
        return self.model.predict(x2)

    @property
    def feature_importances_(self):  # 特征重要性（增益权重）
        if self.model is None:
            raise RuntimeError("The model must be trained with fit () method first")
        return self.model.get_booster().get_score(importance_type='weight')


"""
(1)在梯度提升树（GBDT）（例如 LightGBM）中，默认情况下，每棵树是根据前一棵树的残差来训练的，而训练数据集是固定的，不使用自助采样（bootstrap）。
因此，每棵树/弱学习器（n_estimators）上的样本不会重复。在训练每棵树时，所有的样本都会参与训练，但每棵树的训练数据集是通过拟合前一棵树的残差来构建的。
以LightGBM为例，max_depth控制树的最大深度为6时，连续分裂6次，最多可有2^6-1=63个节点数，2^6=64个理论最大叶子数num_leaves，实际中会更少，
例如min_data_in_leaf 设置了每个叶子节点中最少需要的样本数为5，190个样本则最多有38个节点，当深度为5的时候，节点数已经到31了，
所以num_leaves和max_depth要相互配合。
(2)在集成方法中，使用了bootstrap方法（即自助采样法）来生成每棵树的训练数据。
在这种方法中，每棵树的训练数据是通过从原始数据集中随机抽样来构建的，而抽样是有放回的，意味着一个样本可以被重复选中，同时也有可能某些样本不会被选中。
维度	        XGBoost exact	LightGBM
计算精度	    高（精确分割）	    中（直方图近似）
内存消耗	    高	            低（约XGBoost的1/6）
大数据适应性	差	            优秀
逼近方案	    -	            设置max_bin=65535
"""

# 构建LightGBM回归模型+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
lgbm = EarlyStoppingLGBMRegressor(
    objective='regression',  # 回归任务，multiclass表示完成多分类任务
    # num_leaves=6,  # 每棵树的最大叶子节点数，控制模型复杂度,要小于2^max_depth,这组数据比较小，不需要太多叶子节点 [4,7]
    # learning_rate=0.1,  # 默认0.1  [0.06 0.12]
    n_estimators=150,  # 迭代的次数（树的数量），越大训练效果越好，但会增加计算时长，不过LightGBM已经极快了
    # max_bin=20,  # 默认为255，连续特征在训练时会被离散化为多个bin，更多的bin可以捕捉更多的特征细节，也会增加计算复杂度,[15, 30]
    # bagging_fraction=1.0,  # 样本抽样：每棵树训练时使用的数据样本比例，小于（默认值）1.0时开启随机抽样机制，比例小则计算效率高
    # bagging_freq=10,  # 每10次迭代进行一次数据样本的随机抽样，默认为0（不启用），要和bagging_fraction配合使用
    # bagging_seed=42,  # 随机种子
    max_depth=6,  # 树的深度,会限制num_leaves的数量
    # lambda_l1=0.1,  # L1 正则化, 惩罚复杂模型(将一些特征的权重设为零，防止过拟合)
    lambda_l2=0.1,  # L2 正则化（使部分权重变得较小）
    # feature_fraction=0.8,  # 特征抽样：每棵树训练时使用的特征比例。默认1.0，全部参与
    min_data_in_leaf=10,  # 每个叶子节点至少需要 10 个样本
    # verbosity=1,  # 输出训练日志
    importance_type='gain',  # 信息增益（Gain）,在决策树的分裂过程中，选择某个特征进行分裂时，损失函数（如均方误差、交叉熵）的减少量。gain 值越大，表明该特征对模型预测的贡献越大
    # 若需直接对比特征的总贡献，使用 gain; 若需分析特征被使用的频率，使用 split。
    early_stopping_rounds=6,  # 这个和下一行的val_size是LGBMRegressor()中不自带的
    val_size=0.1,
    random_state=123  # 可以确保每次训练时，模型的初始化和其他随机过程具有相同的起始状态，从而保证模型训练的结果可重复
)
param_grid1 = {  # 修改参数名前缀为 model__，确保与 Pipeline 结构一致。regressor是为了配合反归一化
    'model__regressor__max_bin': [20, 25, 30, 35],  # 尝试不同的 max_bin 设置
    'model__regressor__num_leaves': [2, 4, 6, 8],  # 其他超参数可以一起调优
    'model__regressor__learning_rate': [0.1, 0.12, 0.15, 0.2],
}
# combination1 = Pipeline([('scaler', MinMaxScaler()), ('model', lgbm)])  # 仅仅进行特征归一化，而不进行标签归一化
# TransformedTargetRegressor在训练阶段对目标变量应用一个指定的变换函数（func），在预测阶段使用相应的逆变换函数（inverse_func）恢复原始尺度的目标值。
# 这样可以在回归模型内部处理目标变量的非线性关系或尺度问题，而无需手动对数据进行变换。
combination1 = Pipeline([
    ('scaler_X', MinMaxScaler()),  # 归一化特征
    ('model', TransformedTargetRegressor(
        regressor=lgbm,
        transformer=MinMaxScaler()  # 归一化标签
    ))
])
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 构建XGBoost回归模型++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
xgb_model = XGBoostEarlyStoppingRegressor(  # 显式传入的参数优先级最高，会覆盖类内部定义的默认值或通过**kwargs传递的参数。
    # max_depth=2,
    # learning_rate=0.1,  # 默认0.3 ,越小速度越慢,稳定性越强
    n_estimators=150,  # 迭代次数（树的数量），跟LightGBM给一样的
    objective='reg:squarederror',  # 默认均方误差MSE为回归目标函数，训练模型时最小化损失，二分类：binary:logistic，多分类：multi:softmax
    eval_metric='rmse',  # 评估指标，和reg:squarederror不冲突
    booster='gbtree',  # 使用基于决策树的梯度提升模型
    # gamma=0,  # 分裂时最小损失减少值，值越大分裂越保守，树更简单，默认为0
    # min_child_weight=1,  # 每个叶子节点的最小权重,默认为1，值越大越保守，防止过拟合；对不平衡数据可适当增大
    # subsample=0.9,  # 每棵树训练所用的样本比例
    # colsample_bytree=0.8,  # 每棵树使用的特征比例
    reg_lambda=0.1,  # L2 正则化
    # reg_alpha=0,  # L1 正则化
    importance_type='gain',  # 3个特征重要性：'weight'（默认）：基于特征被用作分裂节点的次数；'gain'：基于特征带来的平均信息增益；'cover'：基于特征覆盖的样本数。
    # 若需对比特征的平均贡献，使用 gain; 若需分析特征被使用的频率，使用 weight。
    early_stopping_rounds=6,  # 早停，避免过拟合
    tree_method='exact',  # hist（直方图近似）加速训练，适合大数据；exact 精确贪心适合小数据
    random_state=123,  # 设置随机种子，模型的初始化和其他随机过程具有相同的起始状态
)
param_grid2 = {
    'model__regressor__max_depth': [2, 4, 6, 8],
    'model__regressor__learning_rate': [0.1, 0.15, 0.2, 0.25],
}
combination2 = Pipeline([
    ('scaler_X', MinMaxScaler()),  # 归一化特征
    ('model', TransformedTargetRegressor(
        regressor=xgb_model,
        transformer=MinMaxScaler()  # 归一化标签
    ))
])
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 构建随机森林模型++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
forest = RandomForestRegressor(
    # n_estimators=200,  # 树（弱学习器）的数量,不太好自己实现早停（原生基于并行计算），所以需要交叉验证
    # max_depth=10,  # 树的最大深度，默认为None(无深度限制，一直生长到不会分裂（拟合精度增加），仅对RF适用)
    random_state=42,  # 控制固定随机数，确保模型可复现（对于大规模数据集，不设置的话精度变化较小）
    min_samples_split=2,  # 内部节点再划分所需的最小样本数，默认为2
    min_samples_leaf=2,  # 最小叶子节点样本数，默认为1（每个叶子节点至少包含一个样本，防止过拟合）
    bootstrap=True,  # 是否使用自助法（bootstrap）,默认为True
    n_jobs=1  # 用于并行计算的CPU核心数，-1表示所有可用的CPU核心,交叉验证的时候全用上，后面往框架里套的时候多种方案用相同的核心数
)
# 基尼不纯度（Gini impurity）：对于分类任务，特征重要性通常基于基尼不纯度的减少。 均方误差（MSE）：对于回归任务，特征重要性通常基于均方误差的减少。
param_grid3 = {
    'model__regressor__max_depth': [4, 6, 8, 10],
    # 'model__regressor__min_samples_split': [2, 3, 5],
    'model__regressor__n_estimators': [50, 75, 100, 125],
    # 'model__regressor__min_samples_leaf': [1, 2, 3],
}
combination3 = Pipeline([
    ('scaler_X', MinMaxScaler()),  # 归一化特征
    ('model', TransformedTargetRegressor(
        regressor=forest,
        transformer=MinMaxScaler()  # 归一化标签
    ))
])


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 构建MLP模型，动态构建模型的函数+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 对权重（kernel）施加 L1 和 L2 正则化，通过加入惩罚项来防止过拟合。l1=0.0001 表示 L1 正则化项的系数，l2=0.001 表示 L2 正则化项的系数。
# bias_regularizer=l2(0.005)：对偏置项进行 L2 正则化，以减少模型的过拟合。
def build_ann_model(hidden_units=16):
    model = models.Sequential([
        layers.Input(shape=(13,)),  # 输入层，特征数
        layers.Dense(hidden_units, activation='relu',  # 激活函数为ReLU 并不是越大越好，太大可能会过拟合，泛化能力下降
                     kernel_regularizer=l1_l2(l1=0.0001, l2=0.001), bias_regularizer=l2(0.005)),
        Dropout(0.2, noise_shape=None, seed=None),  # 接在激活层后, 这一层随机丢弃 10% 的神经元，以减少模型的过拟合
        # layers.Dense(hidden_units // 2, activation='relu'),  # 第二层为第一层的一半
        layers.Dense(hidden_units, activation='relu'),  # 第二层为第一层的一半
        Dropout(0.2 * 0.8, noise_shape=None, seed=None),  # 接在激活层后
        layers.Dense(1),  # 输出层，单个神经元，线性激活函数（默认）
    ])
    model.compile(
        optimizer=optimizers.Adam(learning_rate=0.01),  # Adam优化器学习率默认为0.001，越大收敛越快，但同时不稳定，误差会大一些
        loss='mean_squared_error'  # 损失函数MSE
    )
    return model


# 定义学习率调整策略
reduce_lr = ReduceLROnPlateau(
    monitor='val_loss',  # 监控验证集损失,较低的 val_loss值说明模型在验证集上预测的结果与真实标签之间的差异较小，模型表现较好,MAE也可以作为监控指标
    factor=0.5,  # 学习率减少到原来的 50%
    patience=3,  # 在 3 个 epoch 内没有改善时调整学习率
    min_lr=1e-5,  # 最低学习率 ,第三个早停机制
    verbose=0  # 输出调整学习率的信息
)
# 定义早停策略
early_stopping_ANN = EarlyStopping(
    monitor='val_loss',  # 监控验证集损失
    patience=6,  # 连续6个epoch精度不改善时中断训练
    min_delta=0.001,  # 另一个早停机制，当验证损失改善小于 0.001 时，训练将停止
    verbose=0,  # 输出早停的信息
    restore_best_weights=True  # True 训练停止后，模型将恢复到最佳（即验证集损失最小）时的权重
)
param_grid4 = {
    'model__regressor__model__hidden_units': [8, 16, 32, 64]  # build_ann_model是自建的框架，所以要额外加一个model
}
# 创建Keras回归器包装器
MLP_regressor = KerasRegressor(
    model=build_ann_model,
    epochs=100,  # epochs：遍历模型100次，默认为1
    batch_size=16,  # 并行训练，batch_size默认为32
    validation_split=0.1,
    callbacks=[reduce_lr, early_stopping_ANN],  # 直接传递回调
    verbose=0
)
combination4 = Pipeline([
    ('scaler_X', MinMaxScaler()),  # 归一化特征
    ('model', TransformedTargetRegressor(
        regressor=MLP_regressor,
        transformer=MinMaxScaler()  # 归一化标签
    ))
])
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# 构建k-近邻算法模型++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
knn_model = KNeighborsRegressor(
    # n_neighbors=7,  # 参与预测的最近邻样本数量，默认为5
    # weights='distance',  # 默认为'uniform'，等权加权，'distance'为按距离倒数加权
    algorithm='kd_tree',  # - 'brute'：暴力搜索，'kd_tree'：KD树优化，'ball_tree'：球树优化，‘auto’则会在上面三种算法中做权衡，选择一个拟合最好的最优算法。
    leaf_size=30,  # 影响树结构构建效率的节点容量,默认30，大数据集可以增大
    metric='minkowski',  # 默认闵可夫斯基距离 “minkowski”
    # p=1  # 距离度量参数：1为曼哈顿距离，2为欧氏距离
)
# 基尼不纯度（Gini impurity）：对于分类任务，特征重要性通常基于基尼不纯度的减少。 均方误差（MSE）：对于回归任务，特征重要性通常基于均方误差的减少。
param_grid5 = {
    'model__regressor__n_neighbors': [3, 5, 7, 9],
    'model__regressor__weights': ['uniform', 'distance'],
    'model__regressor__p': [1, 2]
}
combination5 = Pipeline([
    ('scaler_X', MinMaxScaler()),  # 归一化特征
    ('model', TransformedTargetRegressor(
        regressor=knn_model,
        transformer=MinMaxScaler()  # 归一化标签
    ))
])
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 基础数据路径
for idoy in range(181, 182):
    if idoy == 176 or idoy == 182:
        continue
    strdoy = "{:03}".format(idoy)
    path_grid = "B:\\Wide_area_RTPPP\\ML6\\input_data\\features\\data\\" + "input_2024" + strdoy + "_grid.txt"
    path_model = "B:\\Wide_area_RTPPP\\ML6\\input_data\\features\\data\\" + "input_2024" + strdoy + "_model.txt"
    griddata = codecs.open(path_grid, mode='r', encoding='utf-8')
    line0 = griddata.readline()  # 以行的形式进行读取文件
    f2 = codecs.open(path_model, mode='r', encoding='utf-8')
    line2 = f2.readline()

    # 读入文件
    flag, flag2 = 0, 0  # 剔出了一个站
    grid_information, model_information = np.zeros((288 * 382, 16), float), np.zeros((288 * 192, 17), float)
    while line0:
        a = line0.split()
        grid_information[flag, 0:16] = tuple((a[0:16]))  # 将list转换为元组
        line0 = griddata.readline()
        flag = flag + 1
    griddata.close()

    while line2:
        a2 = line2.split()
        model_information[flag2, 0:17] = tuple((a2[0:17]))  # 将list转换为元组
        line2 = f2.readline()
        flag2 = flag2 + 1
    f2.close()

    # EGZTD, VMF3, 5个机器学习模型的ZTD精度，三种决策树的最佳模型的13 * 3个特征重要性, 5个机器学习模型的建模时长,预报时长
    reg_rvmf_rml_fi13 = np.zeros((288, 1 + 1 + 5 + 13 + 13 + 13 + 5 + 5), dtype=float)

    for iepoch in range(110, 288):
        print(iepoch)
        dataall = model_information[iepoch * 192:iepoch * 192 + 192, 2:]

        # 划分训练、测试的特征及标签，13输入 BLH+CMI3+dVMF+Pangu6   1输出dEGZTD
        # input0 = np.hstack((model_information[iepoch * 192:iepoch * 192 + 192, 2:5], model_information[iepoch *
        # 192:iepoch * 192 + 192, 8:9]))
        input0 = np.hstack((model_information[iepoch * 192:iepoch * 192 + 192, 2:8],
                            model_information[iepoch * 192:iepoch * 192 + 192, 8:9] -
                            model_information[iepoch * 192:iepoch * 192 + 192, 15:16],
                            model_information[iepoch * 192:iepoch * 192 + 192, 9:15]))
        label0 = np.hstack((model_information[iepoch * 192:iepoch * 192 + 192, 16:17] -
                            model_information[iepoch * 192:iepoch * 192 + 192, 15:16],
                            model_information[iepoch * 192:iepoch * 192 + 192, 15:16],
                            model_information[iepoch * 192:iepoch * 192 + 192, 8:9]))  # dEGZTD, EGZTD, VMF3
        # label0 = model_information[iepoch * 192:iepoch * 192 + 192, 15:16]
        mask_zero = np.any(dataall == 0.0, axis=1)  # 前两个时间信息不算
        input1 = input0[~mask_zero]  # 在一开始选站的时候已经确定每天的测站缺失的数量不超过10%
        label1 = label0[~mask_zero]  # 除去缺失GNSS-ZTD的历元

        if len(input1) > 170:  # 至少170个训练站
            # 定义5个不同的随机种子
            # random_seeds = [42, 123, 456, 789, 101112]

            # 生成每折的索引分割器（保存训练集和验证集索引）
            custom_cv = []
            flag5 = 0
            rvmf, regztd = np.zeros(5, dtype=float), np.zeros(5, dtype=float)
            kf = KFold(n_splits=5, shuffle=True, random_state=42)
            # 只取第一个分折（确保每种子对应1个分折）
            for train_idx, label_idx in kf.split(input1):
                custom_cv.append((train_idx, label_idx))
                nums = len(label1[label_idx, 2])
                rvmf[flag5] = np.sqrt(
                    np.sum(np.square(
                        label1[label_idx, 2:3] - (label1[label_idx, 0:1] + label1[label_idx, 1:2]))) / nums)  # VMF
                regztd[flag5] = np.sqrt(np.sum(np.square(label1[label_idx, 0:1])) / nums)  # EGZTD
                flag5 += 1
            m_rvmf = np.mean(rvmf)
            m_regztd = np.mean(regztd)
            reg_rvmf_rml_fi13[iepoch, 0] = m_regztd  # 切记最后只统计1小时之后的精度
            reg_rvmf_rml_fi13[iepoch, 1] = m_rvmf

            # break  # KFold.split() 会产生多个训练集和验证集的划分，只关心每个种子下的第一个划分
            # custom_cv = custom_cv[:5]  # 保留前5个划分，一共就5个
            # n_jobs=-1：使用全部可用CPU核心，慎用，不然会崩，n_jobs=1	单进程运行（禁用并行，调试代码），n_jobs=N，指定使用 N 个核心（正整数）；cv=5即为5折交叉验证
            # n_jobs的增加可有效提升实际的训练时间，但mean_fit/score_time等数值会因此增加一点，因为由于 CPU 资源竞争（尤其是计算密集型任务），单个模型训练的墙钟时间（wall time）可能被拉长
            grid_search1 = GridSearchCV(estimator=combination1, param_grid=param_grid1,  # neg_mean_squared_error为负的均方误差
                                        # scoring='neg_mean_squared_error',  # MSE
                                        scoring='neg_root_mean_squared_error',  # 数值越大分越高，相应的，取的绝对值（RMSE）就越小,精度就越高
                                        # scoring='r2',  # 衡量模型对数据的拟合程度
                                        # scoring='neg_mean_poisson_deviance',  # 负泊松偏差
                                        cv=custom_cv, verbose=1, n_jobs=1)
            grid_search2 = GridSearchCV(estimator=combination2, param_grid=param_grid2,  # n_jobs=-1为调用所有CPU，要慎用
                                        scoring='neg_root_mean_squared_error', cv=custom_cv, verbose=1, n_jobs=1)  # 2
            grid_search3 = GridSearchCV(estimator=combination3, param_grid=param_grid3,
                                        scoring='neg_root_mean_squared_error', cv=custom_cv, verbose=1, n_jobs=1)  # 2
            grid_search4 = GridSearchCV(estimator=combination4, param_grid=param_grid4,
                                        scoring='neg_root_mean_squared_error', cv=custom_cv, verbose=1, n_jobs=1)  # 10
            grid_search5 = GridSearchCV(estimator=combination5, param_grid=param_grid5,
                                        scoring='neg_root_mean_squared_error', cv=custom_cv, verbose=1, n_jobs=1)  # 10
            """
            # 对自定义模型进行序列化测试, 确保传递给 GridSearchCV 的数据均为可序列化格式：
            import pickle
            try:
                pickle.dumps(combination1)
            except Exception as e:
                print(f"模型不可序列化: {e}")
            """

            """
            # MinMaxScaler 会按列进行归一化，确保每个特征都在相同的范围内。 这一小部分为不采用交叉验证的调试
            train_index, label_index = list(kf.split(input1))[0]
            input_train, input_test = input1[train_index], input1[label_index]
            out_train, out_test = label1[train_index], label1[label_index]
            scaler = MinMaxScaler()  # 默认缩放到[0 1]之间，其他范围：scaler = MinMaxScaler(feature_range=(-1, 1))
            input_train2 = scaler.fit_transform(input_train)
            input_test2 = scaler.transform(input_test)  # 采用与训练集相同的归一化准则
            scaler2 = MinMaxScaler()
            out_train2 = scaler2.fit_transform(out_train[:, 0:1])  # 训练的时候只用dEGZTD
            out_test2 = scaler2.transform(out_test[:, 0:1])
            forest = RandomForestRegressor(
                n_estimators=200,  # 树（弱学习器）的数量
                max_depth=10,  # 树的最大深度，默认为None(无深度限制，一直生长到不会分裂（拟合精度增加），仅对RF适用)
                random_state=42,  # 控制固定随机数，确保模型可复现（对于大规模数据集，不设置的话精度变化较小）
                min_samples_split=2,  # 内部节点再划分所需的最小样本数，默认为2
                min_samples_leaf=2,  # 最小叶子节点样本数，默认为1（每个叶子节点至少包含一个样本，防止过拟合）
                bootstrap=True,  # 是否使用自助法（bootstrap）,默认为True
                n_jobs=1)  # 用于并行计算的CPU核心数，-1表示所有可用的CPU核心,交叉验证的时候全用上，后面往框架里套的时候多种方案用相同的核心数
            # max_depth=None
            train3_s = time.perf_counter()
            forest.fit(input_train2, out_train2)
            dtrain3 = time.perf_counter() - train3_s
            sim3 = np.transpose(forest.predict(input_test2))  # 预测
            # cim4 = np.transpose(forest.predict(input_train2))
            sim33 = sim3.reshape(-1, 1)  # 重塑为二维数组(inverse_transform的期望数组为二维)：-1表示numpy自动计算行数，1表示只有一个特征
            predict3 = scaler2.inverse_transform(sim33)  # 转置后反归一化
            orms3 = math.sqrt(np.sum((out_test[:, 0:1] ** 2)) / len(predict3))
            rms3 = math.sqrt(np.sum(((predict3 - out_test[:, 0:1]) ** 2)) / len(predict3))
            """
            # time.time()为当前时刻的时间，返回的是自 1970年1月1日以来的秒数（Unix时间戳），time.perf_counter()会包含sleep()休眠时间，适用测量短持续时间
            # train1_s = time.perf_counter()
            # dtrain1 = time.perf_counter() - train1_s  # 结果会偏大一点
            # Scikit-learn 的评分函数遵循以下规则：（1）最大化优化：所有评分函数默认假设值越大越好
            # （2）损失函数转换：对于越小越好的指标（如均方误差 MSE），Scikit-learn 返回其负值（即 neg_mean_squared_error），以便统一通过最大化得分来优化模型
            # 自定义模型类如果未直接暴露LightGBM原生模型的属性，则需通过.model 间接访问,regressor_是对应前面的标签归一化要加入的
            # 为了提供更加可信的feature_importances_结果，模型训练一定要借助交叉验证或者超参数优化器共同执行、以训练一个泛化能力更强的模型
            # feature_importance1 = (best1.named_steps['model']).model.feature_importances_
            grid_search1.fit(input1, label1[:, 0:1])  # 3 * 4 * 3 * 5 = 180； LightGBM  36组超参数验证
            grid_search2.fit(input1, label1[:, 0:1])  # 4 * 4 * 5 = 80； XGBoost 16组超参数验证
            # results1 = grid_search1.cv_results_  # 拟合的信息
            # results2 = grid_search2.cv_results_
            # best1 = grid_search1.best_estimator_  # 最佳模型，以负RMSE为判断指标
            # best2 = grid_search2.best_estimator_

            grid_search3.fit(input1, label1[:, 0:1])  # 4 * 4 * 5 = 80；  Random_forest 16组超参数验证
            grid_search4.fit(input1, label1[:, 0:1])  # 4 * 5 = 20；  ANN/MLP 4组超参数验证
            grid_search5.fit(input1, label1[:, 0:1])  # 4 * 2 * 2 * 5 = 80；  KNN 16组超参数验证

            results1 = grid_search1.cv_results_  # 拟合的信息
            results2 = grid_search2.cv_results_
            results3 = grid_search3.cv_results_
            results4 = grid_search4.cv_results_
            results5 = grid_search5.cv_results_
            best1 = grid_search1.best_estimator_  # 最佳模型，以负RMSE为判断指标
            best2 = grid_search2.best_estimator_
            best3 = grid_search3.best_estimator_
            # best4 = grid_search4.best_estimator_
            # best5 = grid_search5.best_estimator_

            # 计算决策树模型的特征重要性
            feature_i1 = (best1.named_steps['model']).regressor_.feature_importances_  # 特征重要性
            feature_i02 = (best2.named_steps['model']).regressor_.feature_importances_
            feature_i3 = (best3.named_steps['model']).regressor_.feature_importances_
            feature_names = [f'f{i}' for i in range(13)]  # 构建XGboost的特征字典
            # 补全feature_i02中缺失的特征为0
            completed_feature_i2 = {feature: feature_i02.get(feature, 0.0) for feature in feature_names}
            feature_i2 = np.array(list(completed_feature_i2.values()))
            feature_importance1 = feature_i1 / np.sum(feature_i1)  # 转换为百分位数
            feature_importance2 = feature_i2 / np.sum(feature_i2)
            feature_importance3 = feature_i3 / np.sum(feature_i3)  # 即使特征重要性是归一化之后（已经转换为百分位数了）的，这部分也不错，相当于全除以1了

            reg_rvmf_rml_fi13[iepoch, 7:(7 + 13)] = feature_importance1
            reg_rvmf_rml_fi13[iepoch, (7 + 13):(7 + 26)] = feature_importance2
            reg_rvmf_rml_fi13[iepoch, (7 + 26):(7 + 39)] = feature_importance3
            # reg_rvmf_rml_fi13[iepoch, 7 + 39 + 0] = np.sum(results1['mean_fit_time']) * 5  # 180 总训练/建模次数
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 0] = np.mean(results1['mean_fit_time'])  # 180 总训练/建模次数
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 1] = np.mean(results2['mean_fit_time'])  # 80
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 2] = np.mean(results3['mean_fit_time'])  # 80
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 3] = np.mean(results4['mean_fit_time'])  # 20
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 4] = np.mean(results5['mean_fit_time'])  # 80
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 5] = np.mean(results1['mean_score_time'])  # 180 总预测次数
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 6] = np.mean(results2['mean_score_time'])  # 80
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 7] = np.mean(results3['mean_score_time'])  # 80
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 8] = np.mean(results4['mean_score_time'])  # 20
            reg_rvmf_rml_fi13[iepoch, 7 + 39 + 9] = np.mean(results5['mean_score_time'])  # 80
            reg_rvmf_rml_fi13[iepoch, 2] = np.min(np.abs(results1['mean_test_score']))  # 这里的分数就是前面给的负的均方根误差
            reg_rvmf_rml_fi13[iepoch, 3] = np.min(np.abs(results2['mean_test_score']))
            reg_rvmf_rml_fi13[iepoch, 4] = np.min(np.abs(results3['mean_test_score']))
            reg_rvmf_rml_fi13[iepoch, 5] = np.min(np.abs(results4['mean_test_score']))
            reg_rvmf_rml_fi13[iepoch, 6] = np.min(np.abs(results5['mean_test_score']))

    # order = np.linspace(start=0, stop=287, num=288)  # num为分成的个数
    order2 = np.arange(0, 288, 1)  # 第三个为步长
    reg_rvmf_rml_fi13_tra = np.transpose(np.vstack((np.transpose(reg_rvmf_rml_fi13), order2)))
    mask_zero5 = np.all(reg_rvmf_rml_fi13_tra[:, 0:56] == 0.0, axis=1)
    filtered_reg_rvmf_rml_fi13 = reg_rvmf_rml_fi13_tra[~mask_zero5]  # 保留非0行，即为所有格网/测站均有内插值的历元
    np.savetxt("B:\\GOES-16\\result\\CMI_model_" + strdoy + ".txt", filtered_reg_rvmf_rml_fi13, fmt='%.6f')
