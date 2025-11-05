#!/usr/bin/env python
# coding: utf-8

# 缺少的库直接安装pip/conda
import math
import codecs
import warnings
import numpy as np
import pandas as pd
import sklearn.metrics as metrics
import seaborn as sns
# import xgboost as xgb  # 原生API，对应xgb_model = xgb.XGBRegressor()
from xgboost import XGBRegressor
from geopy.distance import geodesic, great_circle  # 计算球面距离
from sklearn.ensemble import RandomForestRegressor
from sklearn.base import BaseEstimator, RegressorMixin  # 这两个类是scikit-learn中的基类，确保自定义的模型类符合scikit-learn的接口标准
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler, PolynomialFeatures
from sklearn.model_selection import train_test_split, KFold, GridSearchCV
from sklearn.compose import TransformedTargetRegressor
from sklearn.metrics import r2_score
from scikeras.wrappers import KerasRegressor
from scipy.special import factorial, gammaln, hyp2f1, lpmv  # 阶乘，伽马函数，超几何函数, 连带勒让德函数
from scipy.optimize import root_scalar, brentq, newton, curve_fit


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

nk_yao = np.array([[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [6.380, 4.840, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [10.49, 10.49, 8.360, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [15.31, 14.79, 14.26, 11.69, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [19.60, 19.60, 18.75, 17.86, 14.93, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [24.29, 23.97, 23.64, 22.53, 21.36, 18.13, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [28.65, 28.65, 28.09, 27.52, 26.20, 24.79, 21.29, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [33.28, 33.04, 32.81, 32.05, 31.28, 29.78, 28.18, 23.73, 0.000, 0.000, 0.000, 0.000, 0.000],
                   [37.67, 37.67, 37.25, 36.83, 35.91, 34.97, 33.30, 31.52, 27.55, 0.000, 0.000, 0.000, 0.000],
                   [42.27, 42.09, 41.90, 41.32, 40.74, 39.68, 38.58, 36.78, 34.83, 30.65, 0.000, 0.000, 0.000],
                   [46.69, 46.69, 46.35, 46.01, 45.30, 44.57, 43.38, 42.15, 40.21, 38.11, 33.74, 0.000, 0.000],
                   [51.27, 51.12, 50.96, 50.49, 50.02, 49.19, 48.33, 47.03, 45.68, 43.61, 41.37, 36.82, 0.000],
                   [55.70, 55.70, 55.41, 55.13, 54.54, 53.95, 53.01, 52.05, 50.64, 49.17, 46.98, 44.62, 39.89]])


def muller_initial_guesses(lll, mmm):
    delta_min = 0.1  # 最小步长阈值
    # 渐近公式初始化
    if lll <= 5:
        lambda_approx = lll * (lll + 1)
        delta = max(0.1 * lll, delta_min)  # lll=0的时候，取delta_min，不然后导致后面muller的初始值全为0
    else:
        lambda_approx = (lll + 0.5) ** 2 - mmm ** 2
        delta = max(0.1 * (lll - mmm), 0.5)

    # 生成三个初始猜测值
    x0 = lambda_approx - delta
    x1 = lambda_approx
    x2 = lambda_approx + delta
    return x0, x1, x2


def muller_method(fun, x0, x1, x2, tol=1e-10, max_iter=1000):  # 马勒求根法，用于辅助非整数勒让德阶数的求解
    """
    Muller方法求根实现
    :param fun: 目标函数
    :param x0, x1, x2: 穆勒法所需的三个初始猜测点。这些点不应完全相等，且它们应该围绕根的区域
    :param tol: 容差表示算法在何种程度上认为求得的根足够准确。默认值为1e-10,可以自己设置，不过这已经够高了
    :param max_iter: 最大迭代次数，防止进入无限循环，默认值为1000，可适当调整,可以多给点，不然不收敛
    :return: 近似根
    """
    # 初始三点检查
    if x0 == x1 or x1 == x2 or x0 == x2:
        raise ValueError("初始点必须互不相同")

    for _ in range(max_iter):  # 添加了一个最大迭代次数（max_iter=100）限制，防止出现无限循环的问题
        f0, f1, f2 = fun(x0), fun(x1), fun(x2)

        # 系数计算部分（严格遵循理论推导）,重新定义相对于x2的偏移量
        h0 = x0 - x2
        h1 = x1 - x2
        delta0 = f0 - f2  # 理论推导中的delta定义,直接使用函数差值
        delta1 = f1 - f2

        # 计算公共分母，确保分母不为0
        denominator = h0 * h1 * (h0 - h1)  # 确保其中任意一个不为0，且二者不相等
        if denominator == 0:
            raise RuntimeError("三点共线导致分母为零，无法继续迭代")

        # 严格按理论公式计算系数
        a = (delta0 * h1 - delta1 * h0) / denominator
        b = (delta1 * h0 ** 2 - delta0 * h1 ** 2) / denominator
        c = f2  # P(x2) = f2

        # 根的计算（修正符号问题）
        discriminant = b ** 2 - 4 * a * c

        # 处理复数根,在球冠谐函数问题中，若球冠半角（如18.5°）与阶数l,m的组合导致方程无实根，复根会成为唯一解
        if discriminant < 0:  # 根号下的负数
            sqrt_d = np.sqrt(-discriminant) * 1j  # 1j代表的是虚数部分为 1 的复数
        else:
            sqrt_d = np.sqrt(discriminant)

        # 构造两个可能的解
        denom_plus = b + sqrt_d
        denom_minus = b - sqrt_d

        # 选择分母绝对值较大的解（数值稳定性）
        if abs(denom_plus) > abs(denom_minus):
            h = -2 * c / denom_plus  # 注意：这里的负号来自理论公式推导
        else:
            h = -2 * c / denom_minus
        x3 = x2 + h  # 修正符号：理论公式要求 x2 + h

        # 收敛性检查
        # 检查迭代步长h（即当前近似解x3与前一点x2的差值）是否小于设定的容差tol。
        # 意义：若步长足够小，说明近似解变化已趋于稳定，可认为收敛到真实根，终止迭代。
        if abs(h) < tol:
            # 返回实数解（如果虚部可忽略）
            if abs(x3.imag) < 1e-15:  # x3复数的虚部的绝对值，如果非常小，认为它几乎为零，意味着x3近似为一个实数
                return x3.real  # 直接返回x3的实部
            return x3

        # 更新策略：抛弃最远的点（标准Muller方法更新方式）,计算三个点到新点的距离
        distances = [abs(x3 - x) for x in [x0, x1, x2]]
        # 找到最远点的索引
        furthest = np.argmax(distances)  # 取出distances中元素最大值所对应的索引

        # 更新三个点（保留两个最近点+新点）,用x3替换最远的
        if furthest == 0:
            x0, x1, x2 = x1, x2, x3
        elif furthest == 1:
            x0, x1, x2 = x0, x2, x3
        else:
            x0, x1, x2 = x0, x1, x3

    raise RuntimeError(f"未在 {max_iter} 次迭代内收敛")


def find_nk(mm, kk, theta0, mode):  # 通过球冠半角求解非整数勒让德阶数
    """
    任意二阶线性常微分方程，都可以化成所谓的Sturm-Liouvile方程，简称S-L型方程，S-L方程中含有参数，在一定的边界条件下，只有取某些特定的值时，
        方程才满足边界条件的非零解，这种值称为问题的本征值，而相应方程的解称为本征函数。
    根据 k-m 的奇偶性选择方程求解非整数阶 n_k(m)，球冠模型半角有关，可由 Dirichlet边界条件、Neumann边界条件和混合边界条件确定
        钱等人，2014；彭等人，2000；测绘学报
    非全局区域对边界条件的额外要求，是Sturm-Liouville理论在受限几何空间中的必然体现，Sturm-Liouville理论的核心要求，其解的正交完备性依赖于算子的自伴性（即Hermitian性质），
        而这需要满足特定边界条件。对球谐函数而言，对应的微分算子是角向拉普拉斯算子的本征问题，其自伴性要求边界处函数值或导数满足特定关系。
        在物理学中, 有一类边界条件是物理系统所处的具体物理环境的数学反映, 它需要人们具体规定, 这种边界条件称为“人为边界条件”, 简称边界条件.
            人为截断点处的条件需反映实际物理约束：
            Dirichlet条件（场值归零）：模拟场在边界处消失（如局部重力异常在区域外衰减）
            Neumann条件（法向导数归零）：模拟通量守恒（如地热流在边界的连续性）
        还有一类边界条件是物理系统本身的固有性质的数学反映, 它是客观存在的, 与物理环境无关, 称其为“自然边界条件”, 其数学表现为要求解函数单值, 有界等.
    （1）全局情况（自然端点θ∈[0,π]） ，如全球VTEC、重力场建模
        在完整球面上，极角θ的端点θ=0和θ=π处，解的自然行为（如有限性）已隐含了边界条件。例如，Legendre函数在θ=0和θ=π处自动满足正则性（无奇点），无需显式施加额外条件即可保证算子自伴性。
    （2）非全局情况（人为截断点，如球冠θ≤θ₀），如区域水汽、ZTD模型
        当区域截断至θ=θ₀时，θ₀成为新的边界点。此时必须显式施加边界条件（如Dirichlet条件Y(θ₀)=0或Neumann条件dY/dθ|θ₀=0），否则算子将不再自伴，导致：
        a.本征值可能为复数 b.本征函数失去正交性 c.解空间不完备，无法展开任意函数
    参数：
        m: 级/次数
        k: 阶数（k >= m）
        theta0: 球冠半角（弧度）
        initial_guess: 初始猜测值

    返回：
        n_k(m) 的根
    """
    t0 = np.cos(theta0)  # 自变量,球冠半角
    is_even = (kk - mm) % 2 == 0  # 判断奇偶性,kk为阶数（主导），mm为次数（伴随）

    def equation(ll):  # 这里ll表示非整数勒让德阶数，为待求解
        # 定义超几何函数参数
        aa = mm - ll
        bb = mm + ll + 1
        cc = mm + 1
        zz = (1 - t0) / 2

        if is_even:  # kk - mm 为偶数，导数在边界θ0处为零
            # Neumann条件
            # 计算超几何函数及其导数,hyp2f1可解超几何函数
            # dF_dz = (aa * bb / cc) * hyp2f1(np.float64(aa + 1), np.float64(bb + 1), np.float64(cc + 1), np.float64(zz))
            # dP_dtheta = np.sin(theta0) * dF_dz / 2  # 链式法则计算导数,这里少了一部分导数
            F_n = hyp2f1(np.float64(aa), np.float64(bb), np.float64(cc), np.float64(zz))
            if np.iscomplex(aa):
                print(f"复数出现：aa={aa}, bb={bb}, cc={cc}, zz={zz}")

            # 计算 F(n_k - 1, m, t0)
            aa_prev = mm - (ll - 1)  # 上一个aa值
            bb_prev = (ll - 1) + mm + 1  # 上一个bb值
            F_n_minus_1 = hyp2f1(np.float64(aa_prev), np.float64(bb_prev), np.float64(cc), np.float64(zz))
            return ll * t0 * F_n - (ll - mm) * F_n_minus_1
        else:  # kk - mm 为奇数
            # Dirichlet条件：函数值在θ0处为0
            F_n = hyp2f1(np.float64(aa), np.float64(bb), np.float64(cc), np.float64(zz))  # 方程变为 F = 0，求解使函数值为零的 n_k(m)
            return F_n

    # 初始猜测点：k-0.1, k, k+0.1
    # xx0, xx1, xx2 = muller_initial_guesses(kk, mm)  # 方法一
    try:
        if mode == 0:
            # x_int = (mm + 0.5) / np.sin(theta0)  # 方法二
            x_int = kk - 0.5 + np.sqrt(np.square(kk + 0.5) - kk * (kk + 1) * np.square(np.sin(theta0)))
            root = muller_method(
                equation, x0=x_int - 0.5, x1=x_int, x2=x_int + 0.5, tol=1e-10, max_iter=1000)  # 在整阶的基础上调整最优值
        elif mode == 1:
            root = muller_method(
                equation, x0=nk_yao[kk, mm] - 0.5, x1=nk_yao[kk, mm], x2=nk_yao[kk, mm] + 0.5, tol=1e-10, max_iter=1000)
        else:
            root = muller_method(
                equation, x0=nk_yao[kk, mm] - 0.1, x1=nk_yao[kk, mm], x2=nk_yao[kk, mm] + 0.1, tol=1e-10, max_iter=1000)
        # root = muller_method(equation, x0=nk_yao[kk, mm] - 0.5, x1=nk_yao[kk, mm], x2=nk_yao[kk, mm] + 0.5,
        #                      tol=1e-10)  # 在整阶的基础上调整最优值
        # root = brentq(equation, kk * np.sin(theta0) - 0.5, kk * np.sin(theta0) + 0.5, xtol=1e-12)
        # root = newton(equation, x0=nk_yao[kk, mm])  # 提供一个初始值
    except RuntimeError as e:
        print(f"警告：{e}, m={mm}, k={kk}")
        root = kk  # 返回默认值
    return root


def compute_associated_legendre(ll, mm, theta0, eps=1e-14, mode=1):  # 通过非整数勒让德阶数、勒让德级数求函数值
    # ll: 实数（小数）阶数（非整数勒让德级数,就是nk(m)）, mm: 整数级数（次数） (m >= 0), theta: 极角（弧度）,eps: 级数截断精度
    cos_theta = np.cos(theta0)  # np.clip将其限制在 [0, 1] 的范围内,应对极小的的负数（例如−1e−16），一般不会出现
    sin_theta = np.sqrt(np.clip(1 - cos_theta ** 2, 0.0, 1.0))  # 处理θ范围,通过余弦值求正弦值
    # 在球坐标系中，极角θ 表示点与极轴（通常为 z 轴）的夹角，范围为[0,π]。 在此范围内，sinθ≥0,因此无需考虑负号。

    # 步骤1：施密特正交化，计算归一化/规格化因子 KK^m_l
    if mm == 0:  # 次数为0
        KK = 1.0
    else:  # 次数大于0
        # 使用对数避免阶乘溢出，例如 100!≈9.3×10^157,超出浮点数范围,不过像是GNSS大气建模（VTEC/ZTD/PWV）中的阶乘数量并不多
        # 伽玛函数（也叫第二类欧拉积分） Γ(𝑧)是对阶乘函数的扩展，对于正整数n  Γ(n)=(n−1)!,gammaln(x) 计算的是伽玛函数的自然对数
        log_num = gammaln(ll + mm + 1) - gammaln(ll - mm + 1)  # 这一部分是中括号里的那部分取对数
        # log_K = (0.5 * np.log(2) - m * np.log(2)
        #          - np.log(factorial(m)) + 0.5 * log_num)  # factorial就是阶乘, 这一句代码和下一部分一样
        log_K = (
                0.5 * np.log(2)
                - mm * np.log(2)
                - gammaln(mm + 1)  # 使用伽马函数gammaln(mm + 1)替代mm!，确保计算稳定性
                + 0.5 * log_num
        )
        KK = np.exp(log_K)  # 取指数，得归一化因子值

    # 下面的while方案可能陷入很大的循环（比如上百万次，因为1e-14很难达到），所以我们设定了另一个循环次数的限制
    # while True:  # 此处ll为非整数勒让德阶数
    #     numerator = (k0 + mm - 1) * (k0 + mm) - ll * (ll + 1)  # 分子
    #     denominator = k0 * (k0 + mm)  # 分母
    #     term_new = (numerator / denominator) * A_ml[k0 - 1]
    #     A_ml.append(term_new)
    #     if abs(term_new) < eps:  # A_ml小于1e-14时中断，这部分term_new会越来越小
    #         break
    #     k0 += 1
    P_ml = 0
    if mode == 0:  # 手动递推公式实现
        # 步骤2：递推计算系数a_k，并累积级数项
        xx = (1 - cos_theta) / 2
        a_k = KK * (sin_theta ** mm)  # 初始项
        xx_k = 1.0  # xx^0
        term = a_k * xx_k
        P_ml = term
        prev_term_abs = float('inf')
        converged = False

        for j0 in range(1, 100):
            # 递推计算a_k，不包含xx的幂次
            numerator = (j0 + mm - 1) * (j0 + mm) - ll * (ll + 1)
            denominator = j0 * (j0 + mm)
            a_k = (numerator / denominator) * a_k  # 更新a_k，不乘xx
            xx_k *= xx  # 更新xx的幂次：xx^k，需要在后续累加的项
            term = a_k * xx_k

            current_term_abs = abs(term)  # 勒让德函数的收敛速度极快，
            P_ml += term

            # 收敛检查
            if current_term_abs < eps:
                converged = True
                break

            # 发散检测：如果当前项不小于前一项，终止循环
            # if current_term_abs >= prev_term_abs:  # 某次迭代之后没变小
            #     warnings.warn(f"发散检测触发于第{k}次迭代。")
            #     break
            # prev_term_abs = current_term_abs
            if not converged:
                warnings.warn(f"未在{max_iter}次迭代内收敛，末项为{current_term_abs:.2e}。")
    if mode == 1:  # 超几何函数方法
        # P_ml = hyp2f1(np.float64(mm - ll), np.float64(ll + mm + 1), np.float64(mm + 1), np.float64((1 - cos_theta) / 2))
        try:
            P_ml = KK * (sin_theta ** mm) * hyp2f1(mm - ll, ll + mm + 1, mm + 1, (1 - cos_theta) / 2)
        except Exception as e:
            print(f"参数错误: a={mm - ll}, b={ll + mm + 1}, c={mm + 1}, z={(1 - cos_theta) / 2}")
            raise

    return P_ml


def pre_compute_nk(K_max, theta0, mode0):  # 这一步是计算所有阶次的n_k(m)，为了给出相应的表格，这个矩阵只由最大阶数及球冠半角决定
    # 预计算所有n_k(m)值，n_matrix是一个下三角矩阵
    # n_0(0)
    # n_1(0) n_1(1)
    # n_2(0) n_2(1) n_2(2) 依次往下类推，这次我算到15阶15次，取最优组合
    n_matrix = np.zeros((K_max + 1, K_max + 1))
    for k in range(K_max + 1):
        for m in range(k + 1):
            n_matrix[k, m] = find_nk(m, k, theta0, mode=mode0)
    return n_matrix


def build_design_matrix(data_points, R, K_max, n_matrix):  # 构建设计矩阵A
    n_data = len(data_points)  # 样本/观测量数
    n_coeffs = (K_max + 1) ** 2  # 待求系数的个数为（阶数+1）的平方
    A = np.zeros((n_data, n_coeffs))  # 初始化矩阵

    for i in range(n_data):
        theta_s_i, lambda_s_i, r_i = data_points[i]  # 极角（余纬度），方位角（经度角）,测站到地心的距离
        col = 0
        for k in range(K_max + 1):  # 阶数
            for m in range(k + 1):  # 次数
                n_k_m = n_matrix[k, m]  # 非整数勒让德阶数
                P = compute_associated_legendre(n_k_m, m, theta_s_i, eps=1e-14, mode=1)  # 计算连带勒让德函数
                radial_term = (R / r_i) ** (n_k_m + 1)  # 归一化坐标距离
                base = R * radial_term * P  # 系数值
                if m == 0:
                    A[i, col] = base  # cos(0λ) = 1
                    col += 1
                else:  # 求设计矩阵
                    A[i, col] = base * np.cos(m * lambda_s_i)
                    A[i, col + 1] = base * np.sin(m * lambda_s_i)
                    col += 2
    return A


def fit_spherical_cap_coefficients(data_points, obs, R=6378137.0, K_max=12, theta0=np.radians(30), mode00=0):
    # 通过最小二乘法拟合球冠谐系数
    """
    :param data_points: 数据点列表，每个元素为 (r, theta_s, lambda_s)
        r: 数据点到地心的距离（米），三维坐标的平方和开根号，用于计算径向衰减项 (R/r)^[nk(m)+1]
        theta_s: 球冠余纬（弧度），范围 [0,π],表示数据点在球冠坐标系中相对于极点的余纬，用于计算关联勒让德函数
        lambda_s: 球冠经度（弧度），范围 [0,2π),表示数据点在球冠坐标系中相对于极点的经度，用于计算三角函数项cos(mλs)和sin(
    :param obs: 观测值数组，形状为 (n_data,)，这里指的是dZTD
    :param R: 地球半径（米）,默认值为 WGS84 椭球长半轴（6378137.0 米）归一化径向距离，确保R/r无量纲化。
    :param K_max: 最大截断阶数, 阶数越高，分辨率越高，但计算量越大，而且有可能过拟合
    :param theta0: 球冠半角（弧度），默认值为 30度转换后的弧度值,后续根据区域自行修改
    :return: g_k_m, h_k_m 系数矩阵，分别对应余弦项和正弦项的幅值

    Args:
        mode00:
        mode00:
    """
    # 预计算 n_k(m)，非整数勒让德阶数
    n_matrix = pre_compute_nk(K_max, theta0, mode0=mode00)

    # 构建设计矩阵，就是通过所有已知值计算的矩阵
    A = build_design_matrix(data_points, R, K_max, n_matrix)

    # 解最小二乘问题，当然也可以采用(A'A)^-1*A'*obs
    # 可参考 https://numpy.org/devdocs/reference/generated/numpy.linalg.lstsq.html#numpy.linalg.lstsq
    # 正常会返回x、resids、rank、s (最小二乘解、残差、矩阵A的秩，a的奇异值)
    x, _, _, _ = np.linalg.lstsq(A, obs, rcond=None)  # rcond用于处理广义逆矩阵计算时时使用的奇异值截断值，一般不用

    # 将系数向量转换为矩阵，这部分可有可无，也可以直接用x做后续的预测运算
    g_k_m = np.zeros((K_max + 1, K_max + 1))
    h_k_m = np.zeros((K_max + 1, K_max + 1))
    col = 0
    for k in range(K_max + 1):
        for m in range(k + 1):
            if m == 0:
                g_k_m[k, m] = x[col]  # z
                col += 1
            else:
                g_k_m[k, m] = x[col]
                h_k_m[k, m] = x[col + 1]
                col += 2
    return g_k_m, h_k_m


def fit_spherical_cap_coefficients_ridge(data_points, obs, R=6378137.0, K_max=12, theta0=np.radians(30), mode00=0,
                                         alpha=0.0):
    # 通过最小二乘法拟合球冠谐系数
    """
    :param data_points: 数据点列表，每个元素为 (r, theta_s, lambda_s)
        r: 数据点到地心的距离（米），三维坐标的平方和开根号，用于计算径向衰减项 (R/r)^[nk(m)+1]
        theta_s: 球冠余纬（弧度），范围 [0,π],表示数据点在球冠坐标系中相对于极点的余纬，用于计算关联勒让德函数
        lambda_s: 球冠经度（弧度），范围 [0,2π),表示数据点在球冠坐标系中相对于极点的经度，用于计算三角函数项cos(mλs)和sin(
    :param obs: 观测值数组，形状为 (n_data,)，这里指的是dZTD
    :param R: 地球半径（米）,默认值为 WGS84 椭球长半轴（6378137.0 米）归一化径向距离，确保R/r无量纲化。
    :param K_max: 最大截断阶数, 阶数越高，分辨率越高，但计算量越大，而且有可能过拟合
    :param theta0: 球冠半角（弧度），默认值为 30度转换后的弧度值,后续根据区域自行修改
    :return: g_k_m, h_k_m 系数矩阵，分别对应余弦项和正弦项的幅值

    Args:
        mode00:
        mode00:
    """
    # 预计算 n_k(m)，非整数勒让德阶数
    n_matrix = pre_compute_nk(K_max, theta0, mode0=mode00)

    # 构建设计矩阵，就是通过所有已知值计算的矩阵
    A = build_design_matrix(data_points, R, K_max, n_matrix)

    # 应用岭回归正则化，alpha的数值要对比确定，一般以对数形式增长，alpha = [0.001, 0.01, 0.1, 1, 10]
    if alpha > 0:  # 正则化参数alpha大于0时才会进行正则化
        n_coeff = A.shape[1]  # 代表设计矩阵中变量的个数，也就是系数的个数
        # reg_matrix = np.sqrt(alpha) * np.eye(n_coeff)  # np.eye为返回一个对角线为1的单位阵
        # A = np.vstack((A, reg_matrix))  # 正则化矩阵被添加到了设计矩阵的底部，它会在最小二乘求解时对系数施加约束
        reg_diag = np.zeros(n_coeff)
        reg_diag[1:] = np.sqrt(alpha)  # 截距项是设计矩阵的第一列，对截距项不做正则化处理，只对除了第一列之外的进行赋值。
        reg_matrix = np.diag(reg_diag)
        A = np.vstack((A, reg_matrix))
        obs = np.vstack((obs, (np.zeros(n_coeff)).reshape(-1, 1)))  # 添加一个零向量来保持矩阵维度的正确性，本质上增加了一些0样本

    # 解最小二乘问题，当然也可以采用(A'A)^-1*A'*obs
    # 可参考 https://numpy.org/devdocs/reference/generated/numpy.linalg.lstsq.html#numpy.linalg.lstsq
    # 正常会返回x、resids、rank、s (最小二乘解、残差、矩阵A的秩，a的奇异值)
    x, _, _, _ = np.linalg.lstsq(A, obs, rcond=None)  # rcond用于处理广义逆矩阵计算时时使用的奇异值截断值，一般不用

    # 将系数向量转换为矩阵，这部分可有可无，也可以直接用x做后续的预测运算
    g_k_m = np.zeros((K_max + 1, K_max + 1))
    h_k_m = np.zeros((K_max + 1, K_max + 1))
    col = 0
    for k in range(K_max + 1):
        for m in range(k + 1):
            if m == 0:
                g_k_m[k, m] = x[col]  # 这个可能就是常数项
                col += 1
            else:
                g_k_m[k, m] = x[col]
                h_k_m[k, m] = x[col + 1]
                col += 2
    return g_k_m, h_k_m


def predict_dZTD(new_points, g_k_m, h_k_m, R=6378137.0, K_max=12, theta0=np.radians(30), mode00=0):
    # 使用拟合的球冠谐系数预测新位置点的dZTD值
    """
    参数：
        new_points: 新位置点列表，每个元素为 (r, theta_s, lambda_s)
        g_k_m: 拟合的余弦项系数矩阵，形状为 (K_max+1, K_max+1)
        h_k_m: 拟合的正弦项系数矩阵，形状为 (K_max+1, K_max+1)
        R: 地球半径（米）
        K_max: 最大截断阶数
        theta0: 球冠半角（弧度）

    返回：
        dZTD_pred: 预测的PWV值数组，形状为 (n_points,)
    """
    # 预计算非整数阶 n_k(m)
    n_matrix = pre_compute_nk(K_max, theta0, mode0=mode00)  # 这里输入的是球冠半角

    # 初始化预测结果
    n_points = len(new_points)
    dZTD_pred = np.zeros(n_points)

    # 遍历所有新位置点
    for i in range(n_points):
        theta_s, lambda_s, r = new_points[i]
        dZTD = 0.0

        # 遍历所有阶数k和次数m
        for k in range(K_max + 1):
            for m in range(k + 1):
                n_k_m = n_matrix[k, m]

                # 计算关联勒让德函数
                P = compute_associated_legendre(n_k_m, m, theta_s, eps=1e-14, mode=1)

                # 计算径向衰减项
                radial_term = (R / r) ** (n_k_m + 1)

                # 计算基函数值
                base = R * radial_term * P

                # 累加余弦项和正弦项
                if m == 0:
                    dZTD += g_k_m[k, m] * base
                else:
                    dZTD += g_k_m[k, m] * base * np.cos(m * lambda_s)
                    dZTD += h_k_m[k, m] * base * np.sin(m * lambda_s)
        dZTD_pred[i] = dZTD
    return dZTD_pred


# poly
def build_poly_matrix(dB, dL, hh, kk, nn):
    scaled_dL = kk * dL  # 弧度
    scaled_dB = kk * dB
    # include_bias=True：生成的特征矩阵中将包含一个常数项列。这个常数项列的值全部为 1
    poly = PolynomialFeatures(degree=nn, include_bias=True, interaction_only=False)
    A_horizontal = poly.fit_transform(np.column_stack([scaled_dB, scaled_dL]))
    A_full = np.hstack([A_horizontal, hh / 1000, np.power((hh / 1000), 2)])  # 高程转换为km
    return A_full


# poly + exp
def cal_poly_exp_coeff(dB, dL, hh, kk, nn, ZTD_obs):
    scaled_dL = kk * dL  # 弧度
    scaled_dB = kk * dB
    # include_bias=True：生成的特征矩阵中将包含一个常数项列。这个常数项列的值全部为 1
    poly2 = PolynomialFeatures(degree=nn, include_bias=True, interaction_only=False)
    A_horizontal = poly2.fit_transform(np.column_stack([scaled_dB, scaled_dL]))
    # ATA = A_horizontal.T @ A_horizontal  # A_horizontal.T表示A_horizontal的转置矩阵,在 Python 中，使用@来执行矩阵的乘法运算
    # ATb = A_horizontal.T @ ZTD_obs  # @与np.dot()都可以用于矩阵乘法或者向量内积，在高维张量运算的时候建议用np.dot()
    # # np.linalg.solve(A, b) 的作用是找到未知向量，满足矩阵方程A⋅x=b，进而求解系数阵x（A_coeffs）
    # A_coeffs = np.linalg.solve(ATA, ATb)  # 求解线性方程组，ATA为方阵
    A_coeffs, _, _, _ = np.linalg.lstsq(A_horizontal, ZTD_obs, rcond=None)
    ZTD0 = A_horizontal @ A_coeffs  # 水平起算值

    # 垂直改正
    B_vertical = np.hstack([hh / 1000, np.power((hh / 1000), 2), np.ones((hh.shape[0], 1))])
    Beita_coeffs, _, _, _ = np.linalg.lstsq(B_vertical, np.log(ZTD_obs), rcond=None)
    return np.vstack([A_coeffs, Beita_coeffs[0:2, :]])  # 垂直加水平系数


def pre_poly_exp_ztd(dB, dL, hh, kk, nn, coeffs):
    scaled_dL = kk * dL  # 弧度
    scaled_dB = kk * dB
    # include_bias=True：生成的特征矩阵中将包含一个常数项列。这个常数项列的值全部为 1
    poly3 = PolynomialFeatures(degree=nn, include_bias=True, interaction_only=False)
    A_horizontal = poly3.fit_transform(np.column_stack([scaled_dB, scaled_dL]))
    B_vertical = np.hstack([hh / 1000, np.power((hh / 1000), 2)])
    # print("A_horizontal shape:", A_horizontal.shape)
    # print("coeffs[:, :-2] shape:", coeffs[:, :-2].shape)
    # print("B_vertical shape:", B_vertical.shape)
    # print("coeffs[:, -2:] shape:", coeffs[:, -2:].shape)
    ZTD_pre = (A_horizontal @ coeffs[:-2, :]) * np.exp((B_vertical @ coeffs[-2:, :]))
    return ZTD_pre  # 预测的ZTD


# 将水平多项式项与垂直指数项的乘积封装为 curve_fit 兼容的函数形式,这部分有待开发，上面的分两步最小二乘的方案效果不好，可能是误差累积较大
# 下面这部分后续用于AI-weather融合建模的时候再开发、测试一下。而且指数主要针对ZTD/ZWD，补偿量不能确保均为正数，其随高度的变化无法确保为指数变化。
def ztd_model(coor, *params):  # 这部分是构造二阶指数对应的多项式的函数  用curve.fit求解函数的系数

    # params: 待优化参数，顺序为 [水平多项式系数 alpha, 垂直系数 beta1, beta2]
    # 提取特征
    dB, dL, hh = coor  # 位置信息
    kk = 1.0  # 假设 kk=1.0，若需动态调整，需将 kk 作为参数传入

    # 提取参数
    poly = PolynomialFeatures(degree=nn, include_bias=True)
    A_horizontal = poly.fit_transform(np.column_stack([kk * dB, kk * dL]))
    n_alpha = A_horizontal.shape[1]  # 水平多项式系数数量
    alpha = params[:n_alpha]
    beta = params[n_alpha:]

    # 水平项
    horizontal_term = A_horizontal @ alpha  # 计算ZWD0

    # 垂直项
    hh_km = hh.reshape(-1, 1) / 1000
    B_vertical = np.hstack([hh_km, hh_km ** 2])
    vertical_term = np.exp(B_vertical @ beta)

    return (horizontal_term * vertical_term).flatten()


# 基础数据路径
path_test = "B:\\Wide_area_RTPPP\\GNSS\\UDUC_VMF\\test.txt"
path_train_validation = "B:\\Wide_area_RTPPP\\GNSS\\UDUC_VMF\\train_validation.txt"
path_grid_0 = "B:\\Wide_area_RTPPP\\GNSS\\UDUC_VMF\\reference_lat_lon_eh_num_382.txt"
test_information = np.loadtxt(path_test, dtype=np.float64, usecols=range(1, 4))
model_information_0 = np.loadtxt(path_train_validation, dtype=np.float64, usecols=range(1, 4))
grid_information_0 = np.loadtxt(path_grid_0, dtype=np.float64, usecols=range(0, 3))
all_coor_0 = np.vstack((test_information, model_information_0, grid_information_0))

# 左纬度、右经度
max_lat = np.max(all_coor_0[:, 0])
min_lat = np.min(all_coor_0[:, 0])
max_lon = np.max(all_coor_0[:, 1])
min_lon = np.min(all_coor_0[:, 1])

box_lat_lon = np.array([[max_lat, max_lon],  # 东北/右上角
                        [max_lat, min_lon],  # 西北/左上角
                        [min_lat, max_lon],  # 东南/右下角
                        [min_lat, min_lon]])  # 西南/左下角
# 球冠极点（球冠中心点）的选取，采用地理中心（几何中心）法，可最小化投影变形：
center_lat = (max_lat + min_lat) * 0.5
center_lon = (max_lon + min_lon) * 0.5

# 计算所有边界点与极点的测地距离,以km为单位，6378137，本质上基于WGS84坐标系
# 地球模型：测地线基于椭球体模型（如 WGS84），认为地球是一个赤道略鼓、两极稍扁的椭球体，而非完美球体,geopy 使用 Vincenty 算法（数值迭代法）求解。
distances_geodesic = [geodesic(box_lat_lon[ii1, :], (center_lat, center_lon)).km for ii1 in range(4)]
# 计算所有边界点与极点的大圆距离，great_circle 内部计算已基于 R = 6371.0 km (完美球体)，直接使用球面三角公式，较为快速
# distances_great_circle = [great_circle(box_lat_lon[ii2, :], (center_lat, center_lon)).km for ii2 in range(4)]
# 计算球冠半角（从球心到球冠边缘的最大角度）
SC_half_angle0 = math.degrees(np.max(distances_geodesic) / 6378.137) + 5  # 弧度转为角度
SC_half_angle = np.radians(SC_half_angle0)

# for idoy in range(36, 37):
for idoy in [66, 202]:
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
    # 三种nk求根方案对应12个阶数：（包括差分和直接的ZTD预测），多项式对应k值
    rr_SCHA = np.zeros((288, 4 + 31 * 4), dtype=float)  # 3个可信度指标（RMSE）+3套方案*31个测站的ZTD
    # rr_Poly = np.zeros((288, 15 * 12), dtype=float)
    # frr_Poly = np.zeros((288, 15 * 12), dtype=float)
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

    """
    # 左纬度、右经度
    max_lat = np.max((np.max(grid_information[:, 2]), np.max(model_information[:, 2])))
    min_lat = np.min((np.min(grid_information[:, 2]), np.min(model_information[:, 2])))
    max_lon = np.max((np.max(grid_information[:, 3]), np.max(model_information[:, 3])))
    min_lon = np.min((np.min(grid_information[:, 3]), np.min(model_information[:, 3])))

    # 获取包含 0 的行的索引
    rows_index = np.where(np.any(grid_information[:, 2:3] == 0, axis=1))[0]

    box_lat_lon = np.array([[max_lat, max_lon],  # 东北/右上角
                            [max_lat, min_lon],  # 西北/左上角
                            [min_lat, max_lon],  # 东南/右下角
                            [min_lat, min_lon]])  # 西南/左下角
    # 球冠极点（球冠中心点）的选取，采用地理中心（几何中心）法，可最小化投影变形：
    center_lat = (max_lat + min_lat) * 0.5
    center_lon = (max_lon + min_lon) * 0.5

    # 计算所有边界点与极点的测地距离,以km为单位，6378137，本质上基于WGS84坐标系
    # 地球模型：测地线基于椭球体模型（如 WGS84），认为地球是一个赤道略鼓、两极稍扁的椭球体，而非完美球体,geopy 使用 Vincenty 算法（数值迭代法）求解。
    distances_geodesic = [geodesic(box_lat_lon[ii1, :], (center_lat, center_lon)).km for ii1 in range(4)]
    # 计算所有边界点与极点的大圆距离，great_circle 内部计算已基于 R = 6371.0 km (完美球体)，直接使用球面三角公式，较为快速
    # distances_great_circle = [great_circle(box_lat_lon[ii2, :], (center_lat, center_lon)).km for ii2 in range(4)]
    # 计算球冠半角（从球心到球冠边缘的最大角度）
    SC_half_angle = math.degrees(np.max(distances_geodesic) / 6378.137) + 5  # 弧度转为角度
    """

    # EGZTD, VMF3, 5个机器学习模型的ZTD精度，三种决策树的最佳模型的13 * 3个特征重要性, 5个机器学习模型的建模时长,预报时长
    reg_rvmf_rml_fi13 = np.zeros((288, 1 + 1 + 5 + 13 + 13 + 13 + 5 + 5), dtype=float)

    for iepoch in range(77, 288):
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
        test_input0 = np.hstack((grid_information[iepoch * 382:iepoch * 382 + 382, 2:8],
                                 grid_information[iepoch * 382:iepoch * 382 + 382, 8:9] -
                                 grid_information[iepoch * 382:iepoch * 382 + 382, 15:16],
                                 grid_information[iepoch * 382:iepoch * 382 + 382, 9:15]))  # 虚拟参考站信息
        # test_input0没有0行，所以不需要剔除。

        # label0 = model_information[iepoch * 192:iepoch * 192 + 192, 15:16]
        mask_zero = np.any(dataall == 0.0, axis=1)  # 前两个时间信息不算, 从input0数组中筛选出不包含0.0的行，后续将这些行赋值给input1
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

            grid_search3 = GridSearchCV(estimator=combination2, param_grid=param_grid2,
                                        scoring='neg_root_mean_squared_error', cv=custom_cv, verbose=1, n_jobs=1)  # 2

            grid_search3.fit(input1, label1[:, 0:1])  # 4 * 4 * 5 = 80；  Random_forest 16组超参数验证
            results3 = grid_search3.cv_results_
            best3 = grid_search3.best_estimator_
            grid_res = best3.predict(test_input0)
            # 这个是预测的格网ZTD，在球谐/多项式等最优方法探索的过程中作为虚拟参考站的真值
            grid_ztd = grid_res + grid_information[iepoch * 382:iepoch * 382 + 382, 15:16]

            # 经纬高，VMF3的格网/站基补偿， 预测的格网补偿/真实站基补偿，预测的格网ZTD
            all_ztd1 = np.hstack((grid_information[iepoch * 382:iepoch * 382 + 382, 2:5],  # 仅仅VMF3
                                  grid_information[iepoch * 382:iepoch * 382 + 382, 8:9] -
                                  grid_information[iepoch * 382:iepoch * 382 + 382, 15:16], grid_res, grid_ztd))
            # 仅仅GNSS,包括经纬高，VMF3-ZTD与EGZTD的残差，GNSS-ZTD与EGZTD的残差，GNSS-ZTD
            all_ztd2 = np.hstack((input1[:, 0:3], input1[:, 6:7], label1[:, 0:1], label1[:, 0:1] + label1[:, 1:2]))
            all_ztd = np.vstack((all_ztd1, all_ztd2))  # VMF3+GNSS

            # 将测站的大地坐标(BB, LL)转换为地心坐标 (φ, λ)
            BB = np.radians(all_ztd[:, 0])
            LL = np.radians(all_ztd[:, 1])  # 0.00669438002为WGS-84椭球第一偏心率e的平方，WGS84坐标系下e=sqrt(a^2 - b^2)/a
            phi = np.arctan((1 - 0.00669438002) * np.tan(BB))  # 地心纬度 φ
            lambda_earth = LL  # 地心经度 λ（与大地经度相同）
            # 测试站
            BB_test = np.radians(test_information[:, 0])
            LL_test = np.radians(test_information[:, 1])  # 0.00669438002为WGS-84椭球第一偏心率e的平方，WGS84坐标系下e=sqrt(a^2 - b^2)/a
            phi_test = np.arctan((1 - 0.00669438002) * np.tan(BB_test))  # 地心纬度 φ
            lambda_earth_test = LL_test  # 地心经度 λ（与大地经度相同）

            # 转换为余纬度 θ = π/2 - φ
            theta = np.pi / 2 - phi
            theta_test = np.pi / 2 - phi_test

            # 计算地心直角坐标
            Nr = 6378137.0 / np.sqrt(1 - 0.00669438002 * np.square(np.sin(BB)))  # 卯酉圈曲率半径，** 运算符和 np.sqrt 同样支持数组的逐元素计算
            Xs = (Nr + all_ztd[:, 2]) * np.cos(BB) * np.cos(LL)
            Ys = (Nr + all_ztd[:, 2]) * np.cos(BB) * np.sin(LL)
            Zs = (Nr * (1 - 0.00669438002) + all_ztd[:, 2]) * np.sin(BB)
            distance_earth_center = np.sqrt(np.square(Xs) + np.square(Ys) + np.square(Zs))

            Nr_test = 6378137.0 / np.sqrt(1 - 0.00669438002 * np.square(np.sin(BB_test)))
            Xs_test = (Nr_test + test_information[:, 2]) * np.cos(BB_test) * np.cos(LL_test)
            Ys_test = (Nr_test + test_information[:, 2]) * np.cos(BB_test) * np.sin(LL_test)
            Zs_test = (Nr_test * (1 - 0.00669438002) + test_information[:, 2]) * np.sin(BB_test)
            distance_earth_center_test = np.sqrt(np.square(Xs_test) + np.square(Ys_test) + np.square(Zs_test))

            # 处理球冠极点的地心坐标
            B_N = np.radians(center_lat)
            L_N = np.radians(center_lon)
            phi_N = np.arctan((1 - 0.00669438002) * np.tan(B_N))  # 极点的地心纬度 φ_N
            theta_N = np.pi / 2 - phi_N  # 极点的余纬度 θ_N
            lambda_N = L_N  # 极点的地心经度 λ_N

            # 计算球冠坐标 θ_c 和 λ_c
            cos_theta_c = np.cos(theta_N) * np.cos(theta) + np.sin(theta_N) * np.sin(theta) * np.cos(
                lambda_N - lambda_earth)
            theta_c = np.arccos(cos_theta_c)  # 所有测站的球冠极角（余纬度）不能超过极点的球冠极角（余纬度）
            cos_theta_c_test = np.cos(theta_N) * np.cos(theta_test) + np.sin(theta_N) * np.sin(theta_test) * np.cos(
                lambda_N - lambda_earth_test)
            theta_c_test = np.arccos(cos_theta_c_test)

            # 利用球面三角学公式计算 sin(λ_c) 和 cos(λ_c)
            sin_lambda_c = (np.sin(theta) * np.sin(lambda_earth - lambda_N)) / np.sin(theta_c)
            cos_lambda_c = (np.sin(theta_N) * np.cos(theta) - np.cos(theta_N) * np.sin(theta) * np.cos(
                lambda_earth - lambda_N)) / np.sin(theta_c)
            sin_lambda_c_test = (np.sin(theta_test) * np.sin(lambda_earth_test - lambda_N)) / np.sin(theta_c_test)
            cos_lambda_c_test = (np.sin(theta_N) * np.cos(theta_test) - np.cos(theta_N) * np.sin(theta_test) * np.cos(
                lambda_earth_test - lambda_N)) / np.sin(theta_c_test)

            # 普通 np.arctan(y/x) 无法区分 (x, y) 在不同象限的情况（例如 (1, 1) 和 (-1, -1) 的 y/x 相同，但角度相差 π）。
            # 使用 arctan2 计算 λ_c（单位：弧度）,np.arctan2(sin_lambda_c, cos_lambda_c) 会根据sin_lambda_c和cos_lambda_c的符号自动确定正确的象限。
            lambda_c = np.arctan2(sin_lambda_c, cos_lambda_c)
            lambda_c_test = np.arctan2(sin_lambda_c_test, cos_lambda_c_test)

            # poly，各个测站与中心点的经纬差
            # BB = np.radians(all_ztd[:, 0])
            # LL = np.radians(all_ztd[:, 1])
            # B_N = np.radians(center_lat)
            # L_N = np.radians(center_lon)
            dB = BB - B_N
            dL = LL - L_N
            # dB = all_ztd[:, 0] - center_lat
            # dL = all_ztd[:, 1] - center_lon

            kf_sh = KFold(n_splits=10, shuffle=True, random_state=123)
            # 只取第一个分折（确保每种子对应1个分折）
            train_idx5, label_idx5 = list(kf_sh.split(all_ztd))[0]
            coor_para = np.column_stack((theta_c, lambda_c, distance_earth_center))
            coor_para_test = np.column_stack((theta_c_test, lambda_c_test, distance_earth_center_test))

            all_ztd_VMF3 = all_ztd[0:382, :]
            coor_VMF3 = coor_para[0:382, :]

            all_ztd_GNSS = all_ztd[382:, :]
            coor_GNSS = coor_para[382:, :]

            train_idx6, label_idx6 = list(kf_sh.split(all_ztd_VMF3))[0]
            train_idx7, label_idx7 = list(kf_sh.split(all_ztd_GNSS))[0]
            # np.savetxt("B:\\Wide_area_RTPPP\\ML6\\All_day_accuracy_ML5\\SCHA_Poly_accuracy\\"
            #            + "fit_sta.txt", all_ztd[train_idx5, 0:3], fmt='%.6f', delimiter=" ")  # 用于拟合的测站
            # np.savetxt("B:\\Wide_area_RTPPP\\ML6\\All_day_accuracy_ML5\\SCHA_Poly_accuracy\\"
            #            + "predict_sta.txt", all_ztd[label_idx5, 0:3], fmt='%.6f', delimiter=" ")  # 用于验证最佳阶数的站

            # rr3 = np.zeros((15, 1), dtype=float)

            #  若引入地理对称性条件（例如仅保留偶数次项），系数数量可减少。如果精度达标的化，后续可以考虑以此方法进一步减少模型的系数，但是估计区域情况可能不太适用
            # 方案2和4的可信度指标会虚高
            g_k_m1, h_k_m1 = fit_spherical_cap_coefficients(coor_para[train_idx5, :],
                                                            all_ztd[train_idx5, 4:5], R=6378137.0,  # GNSS测站+ML格网
                                                            K_max=5, theta0=SC_half_angle, mode00=1)
            g_k_m2, h_k_m2 = fit_spherical_cap_coefficients(coor_VMF3[train_idx6, :],  # 仅采用机器学习的格网点的拟合系数
                                                            all_ztd_VMF3[train_idx6, 4:5], R=6378137.0,
                                                            K_max=5, theta0=SC_half_angle, mode00=1)
            g_k_m3, h_k_m3 = fit_spherical_cap_coefficients(coor_GNSS[train_idx7, :],
                                                            all_ztd_GNSS[train_idx7, 4:5], R=6378137.0,  # GNSS测站
                                                            K_max=5, theta0=SC_half_angle, mode00=1)
            g_k_m4, h_k_m4 = fit_spherical_cap_coefficients(coor_para[train_idx5, :],  # VMF3在格网点+站点的拟合系数
                                                            all_ztd[train_idx5, 3:4], R=6378137.0,
                                                            K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq1 = predict_dZTD(coor_para[label_idx5, :],
                                     g_k_m1, h_k_m1, R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq2 = predict_dZTD(coor_VMF3[label_idx6, :],
                                     g_k_m2, h_k_m2, R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq3 = predict_dZTD(coor_GNSS[label_idx7, :],
                                     g_k_m3, h_k_m3, R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq4 = predict_dZTD(coor_para[label_idx5, :],
                                     g_k_m4, h_k_m4, R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            rr_SCHA[iepoch, 0:1] = np.sqrt(np.mean(np.square(dZTD_lsq1.reshape(-1, 1) - all_ztd[label_idx5, 4:5])))
            rr_SCHA[iepoch, 1:2] = np.sqrt(np.mean(np.square(dZTD_lsq2.reshape(-1, 1) - all_ztd_VMF3[label_idx6, 4:5])))
            rr_SCHA[iepoch, 2:3] = np.sqrt(np.mean(np.square(dZTD_lsq3.reshape(-1, 1) - all_ztd_GNSS[label_idx7, 4:5])))
            rr_SCHA[iepoch, 3:4] = np.sqrt(np.mean(np.square(dZTD_lsq4.reshape(-1, 1) - all_ztd[label_idx5, 4:5])))

            # 上面这部分适用于计算可信度指标，后续将ZTD转换为ZWD时，ZHD对于真值和预测值是一样的，因此ZWD的残差和ZTD残差相同，可信度指标一样
            # 这部分为计算经验模型的补偿量
            g_k_m11, h_k_m11 = fit_spherical_cap_coefficients(coor_para, all_ztd[:, 4:5], R=6378137.0,
                                                              K_max=5, theta0=SC_half_angle, mode00=1)
            g_k_m22, h_k_m22 = fit_spherical_cap_coefficients(coor_VMF3, all_ztd_VMF3[:, 4:5], R=6378137.0,
                                                              K_max=5, theta0=SC_half_angle, mode00=1)
            g_k_m33, h_k_m33 = fit_spherical_cap_coefficients(coor_GNSS, all_ztd_GNSS[:, 4:5], R=6378137.0,
                                                              K_max=5, theta0=SC_half_angle, mode00=1)
            g_k_m44, h_k_m44 = fit_spherical_cap_coefficients(coor_para, all_ztd[:, 3:4], R=6378137.0,
                                                              K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq11 = predict_dZTD(coor_para_test, g_k_m11, h_k_m11,  # fus1
                                      R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq22 = predict_dZTD(coor_para_test, g_k_m22, h_k_m22,  # fus2
                                      R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq33 = predict_dZTD(coor_para_test, g_k_m33, h_k_m33,  # GNSS(only)
                                      R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            dZTD_lsq44 = predict_dZTD(coor_para_test, g_k_m44, h_k_m44,  # VMF3(only)
                                      R=6378137.0, K_max=5, theta0=SC_half_angle, mode00=1)
            for ik in range(0, 31):
                rr_SCHA[iepoch, 4 + ik] = dZTD_lsq11[ik]
                rr_SCHA[iepoch, 4 + ik + 31] = dZTD_lsq22[ik]
                rr_SCHA[iepoch, 4 + ik + 62] = dZTD_lsq33[ik]
                rr_SCHA[iepoch, 4 + ik + 93] = dZTD_lsq44[ik]

            # reg_rvmf_rml_fi13[iepoch, 7:(7 + 13)] = feature_importance1

    # order = np.linspace(start=0, stop=287, num=288)  # num为分成的个数
    order2 = np.arange(0, 288, 1)  # 第三个为步长
    reg_rr_SCHA = np.transpose(np.vstack((np.transpose(rr_SCHA), order2)))
    mask_zero1 = np.all(reg_rr_SCHA[:, 0:72] == 0.0, axis=1)
    filtered_reg_rr_SCHA = reg_rr_SCHA[~mask_zero1]  # 保留非0行，即为所有格网/测站均有内插值的历元
    np.savetxt("B:\\Wide_area_RTPPP\\ML6\\All_day_accuracy_ML5\\SCHA_Poly_accuracy\\scheme3_." + strdoy + ".txt",
               filtered_reg_rr_SCHA, fmt='%.6f', delimiter=" ")
