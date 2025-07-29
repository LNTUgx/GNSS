#!/usr/bin/python

from netCDF4 import Dataset
import matplotlib as mpl

mpl.use('TkAgg')
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from pyproj import CRS, Transformer  # (Coordinate Reference System),pyproj中表示投影的高层次对象，用于定义和处理地理或投影坐标系
from pykrige.ok import OrdinaryKriging  # 克里金的库
from multiprocessing import Pool, Manager
import numpy as np
import codecs
import os


#  https://makersportal.com/blog/2018/11/25/goes-r-satellite-latitude-and-longitude-grid-projection-algorithm#algorithm

# def lat_lon_reproj(nc_folder, nc_indx):
def lat_lon_reproj(nc_folder):
    # os.chdir(nc_folder)
    # full_direc = os.listdir()
    # nc_files = [ii for ii in full_direc if ii.endswith('.nc')]  # 查找以nc为结尾的文件
    # g16_data_file = nc_files[nc_indx]  # select .nc file
    g16_data_file = nc_folder
    print()
    # print(nc_files[nc_indx])  # print file name

    # designate dataset
    g16nc = Dataset(g16_data_file, 'r')
    var_names = [ii for ii in g16nc.variables]
    var_name = var_names[0]
    try:
        band_id = g16nc.variables['band_id'][:]
        band_id = ' (Band: {},'.format(band_id[0])
        band_wavelength = g16nc.variables['band_wavelength']
        band_wavelength_units = band_wavelength.units
        band_wavelength_units = '{})'.format(band_wavelength_units)
        band_wavelength = ' {0:.2f} '.format(band_wavelength[:][0])
        print('Band ID: {}'.format(band_id))  # 输出波段ID
        # print('Band Wavelength: {} {}'.format(band_wavelength, band_wavelength_units))  # 输出波长
    except:
        band_id = ''
        band_wavelength = ''
        band_wavelength_units = ''

    # GOES-R projection info and retrieving relevant constants 都是从GOES-NC文件里读的，可以从Panoply里看到详细内容
    # 后续如果读取GOES-17/18的数据要修改以下等基准
    proj_info = g16nc.variables['goes_imager_projection']  # 官方文件也有相应介绍
    lon_origin = proj_info.longitude_of_projection_origin  # 初始经度 -75.0
    H = proj_info.perspective_point_height + proj_info.semi_major_axis  # perspective_point_height: 35786023.0
    r_eq = proj_info.semi_major_axis  # semi_major_axis: 6378137.0 赤道半径
    r_pol = proj_info.semi_minor_axis  # semi_minor_axis: 6356752.31414  极区半径

    # grid info  “x”和“y”投影坐标仅给出关于每个像素相对于成像仪的角度的信息
    lat_rad_1d = g16nc.variables['x'][:]
    lon_rad_1d = g16nc.variables['y'][:]

    # data info 提取数据变量、单位、影像拍摄截止时间（目前推测为UTC），应用时要采用文件名中显示的播发给用户的时间、以及长文件名
    data = g16nc.variables[var_name][:]
    data_units = g16nc.variables[var_name].units
    data_time_grab = ((g16nc.time_coverage_end).replace('T', ' ')).replace('Z', '')
    data_long_name = g16nc.variables[var_name].long_name

    # close file when finished
    g16nc.close()
    g16nc = None

    # create meshgrid filled with radian angles  只适用于1500*2500的GOES系列
    # 下面的一系列公式是为了. Reprojecting Scan Angles (Radians) to Geographic Coordinates (Lat/Lon)
    lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)

    # lat/lon calc routine from satellite radian angle vectors

    lambda_0 = (lon_origin * np.pi) / 180.0  # 初始经度，角度转弧度
    # a/b/c_var均为过程变量
    a_var = np.power(np.sin(lat_rad), 2.0) + (np.power(np.cos(lat_rad), 2.0) * (
            np.power(np.cos(lon_rad), 2.0) + (((r_eq * r_eq) / (r_pol * r_pol)) * np.power(np.sin(lon_rad), 2.0))))
    b_var = -2.0 * H * np.cos(lat_rad) * np.cos(lon_rad)
    c_var = (H ** 2.0) - (r_eq ** 2.0)

    r_s = (-1.0 * b_var - np.sqrt((b_var ** 2) - (4.0 * a_var * c_var))) / (2.0 * a_var)
    #  我推测(b_var ** 2) - (4.0 * a_var * c_var) 开方后有复数,所以填充为NaN了，导致后面的经纬度及data都在西北角都被mask了(掩码)

    s_x = r_s * np.cos(lat_rad) * np.cos(lon_rad)  # 三个坐标旋转量s_x/y/z可通过像素角度的弦值生成
    s_y = - r_s * np.sin(lat_rad)
    s_z = r_s * np.cos(lat_rad) * np.sin(lon_rad)

    lat = (180.0 / np.pi) * (
        np.arctan(((r_eq * r_eq) / (r_pol * r_pol)) * ((s_z / np.sqrt(((H - s_x) * (H - s_x)) + (s_y * s_y))))))
    lon = (lambda_0 - np.arctan(s_y / (H - s_x))) * (180.0 / np.pi)

    # print test coordinates
    # print('{} N, {} W'.format(lat[318, 1849], abs(lon[318, 1849])))  # 随便输出某个像元的经纬度

    return (lon, lat, data, data_units, data_time_grab, data_long_name,
            band_id, band_wavelength, band_wavelength_units, var_name)
    # 返回了波段及波长等信息


path_grid = "B:\\Wide_area_RTPPP\\GNSS\\UDUC_VMF\\reference_lat_lon_eh_num_undu_382.txt"
path_model = "B:\\Wide_area_RTPPP\\GNSS\\UDUC_VMF\\train_validation.txt"
griddata = codecs.open(path_grid, mode='r', encoding='utf-8')
line0 = griddata.readline()  # 以行的形式进行读取文件
f2 = codecs.open(path_model, mode='r', encoding='utf-8')
line2 = f2.readline()

flag, flag2 = 0, 0
# 经纬高+5个比湿+t2m
blh_num, model_information = np.zeros((382, 2), float), np.zeros((193, 2), float)
while line0:  # 只提取经纬高
    a = line0.split()
    blh_num[flag, 0:2] = tuple((a[0:2]))  # 将list转换为元组
    line0 = griddata.readline()
    flag = flag + 1
griddata.close()

while line2:
    a2 = line2.split()
    model_information[flag2, 0:2] = tuple((a2[1:3]))  # 将list转换为元组
    line2 = f2.readline()
    flag2 = flag2 + 1
f2.close()


def process_file(pathroot, idoy):
    full_direc = os.listdir(pathroot)
    CMI_name = [ii for ii in full_direc if ii.endswith('.nc')]  # 查找nc文件
    # 使用 np.full 初始化一个二维字符串数组，填充值为空字符串
    strings_array = np.full((288, 3), '', dtype='<U90')  # 设置最大字符长度为90

    for irs0 in range(0, len(CMI_name)):
        iband0 = int(((CMI_name[irs0]).split('M6C')[1])[0:2]) - 8  # 8\9\10 波段分别对应0\1\2
        iepoch = (CMI_name[irs0]).split('_c2024')[1]
        hour = int(iepoch[3:5])  # 下面指的是文件名后三位为秒，例如024为2.4s, 596为59.6s，这里除以10后除以60是换算成min
        # 一般一天的数据是从0h4min（传输完毕）开始（文件名里三个时间，分别为开始采集（实时）、结束、传输完毕的时间），按照5min采样，23h59min结束
        # 024这种一中对应的是前一个历元+5min，596是前一个历元+4min（少见的情况），但是转换为min求和再取整正好为5min采样的整数
        min = int(iepoch[5:7]) + round((float(iepoch[7:10]) / 600.0)) + 1  # 这里加1min是为了后面和GNSS数据匹配时前置时间多1min
        strings_array[(hour * 12 + int(min / 5) - 1), iband0] = CMI_name[irs0]
    # 检查每行是否有任何空字符串,~为取反操作，识别非空行.这一步是确保采用的时刻在三个通道都有数据.axis=1：表示按行进行操作
    mask_name = ~np.any(strings_array == '', axis=1)
    # 使用 mask 来选择非空行
    final_array = strings_array[mask_name]
    # 单天数据定义，两个时间信息加三个CMI
    CMI_bt_site = np.zeros((193 * len(final_array), 5), dtype=float)  # 每天每个通道288张影像，5min采样间隔
    CMI_bt_grid = np.zeros((382 * len(final_array), 5), dtype=float)  # 可能会有例如C08有数据，09/10没有数据的情况

    for irs in range(0, len(final_array)):
        print(irs)
        stop, stop2 = 0, 0
        for iband in range(0, 3):
            nc_folder = pathroot + final_array[irs, iband]  # define folder where .nc files are located
            index_band3 = iband  # C08/9/10分别对应0\1\2
            # file_indx = 1  # be sure to pick the correct file. Make sure the file is not too big either,
            # some of the bands create large files (stick to band 7-16)
            iepoch2 = (final_array[irs, iband]).split('_c2024')[1]
            hour2 = int(iepoch2[3:5])
            min2 = int(iepoch2[5:7]) + int(round((float(iepoch2[7:10]) / 600.0)) + 1)

            # main data grab from function above
            lon0, lat0, data0, data_units, data_time_grab, data_long_name, band_id, band_wavelength, band_units, var_name \
                = lat_lon_reproj(nc_folder)

            # np.ma.masked_where根据指定条件（逗号前面的就是条件）将数组中的某些值屏蔽（即掩码）,这一步好像没有必要
            data0 = np.ma.masked_where((data0 == None) | (data0 == '') | np.isnan(np.array(data0, dtype=float)), data0)
            lat0 = np.ma.masked_where((lat0 == None) | (lat0 == '') | np.isnan(np.array(lat0, dtype=float)), lat0)
            lon0 = np.ma.masked_where((lon0 == None) | (lon0 == '') | np.isnan(np.array(lon0, dtype=float)), lon0)
            # data0是MaskedArray类型，在数组中标记某些元素为“屏蔽”（即无效或缺失）（False、true），.filled是填充掩码数组（Masked Array）中掩码元素的函数
            data = data0.filled(9999.9)  # data这个值随便替换，只要能识别就行，pcolormesh输出的时候会掩藏
            lat = lat0.filled(-666.6)  # lat/lon给一个特别小的值，不在图片显示范围内即可
            lon = lon0.filled(-666.6)  # 不能负的太大，不然pcolormesh可能会对极端值进行特殊处理，不能太小，不然会显示在图区内（比如-100）

            geos_proj = CRS.from_dict({  # 这个对应卫星项目自己进行更新
                "proj": "geos",  # 地球静止投影
                "h": 35786023.0,  # 卫星视点高度
                "a": 6378137.0,  # 地球长半轴
                "b": 6356752.31414,  # 地球短半轴
                "lon_0": -75.0,  # 投影原点经度
                "sweep": "x"  # 扫描角轴
            })
            geo_proj = CRS.from_epsg(4326)  # WGS84 坐标系, 创建了一个 CRS 对象，表示 WGS84 坐标系，通常用于全球范围的地理坐标系
            # 因为是地球同步轨道卫星，每个像素相对于成像仪的角度的信息是固定的，因此后续解算的经纬度，xy等都是不变的，每张光学影像的像元的位置固定。

            # 创建转换器  # always_xy=True 确保了传入的经纬度顺序是 (lon, lat)
            transformer = Transformer.from_crs(geo_proj, geos_proj, always_xy=True)
            # 展平经纬度数组进行转换
            flat_lon = lon.ravel()
            flat_lat = lat.ravel()

            # 批量转换经纬度 -> x, y
            flat_x, flat_y = transformer.transform(flat_lon, flat_lat)

            # 将结果重新变为二维数组
            x_coor0 = flat_x.reshape(lon.shape)  # 使用的目标形状与原始数组的形状一致，元素在多维数组中的位置就会恢复到与原始数组对应的位置
            y_coor0 = flat_y.reshape(lat.shape)
            # 将地面测站和格网点的经纬度转换为卫星影像的 x, y 坐标(WGS-84)
            x_grid, y_grid = transformer.transform(blh_num[:, 1], blh_num[:, 0])
            x_site, y_site = transformer.transform(model_information[:, 1], model_information[:, 0])  # 先经度后纬度

            # 这里lat、lon及数据掩膜的可能不是同一个位置，所以要取并集，防止后面选取最近像元的时候，距离满足要求，但是对应位置的亮温数据为空
            valid_points2 = ((~lat0.mask) & (~lon0.mask) & (~data0.mask) & (lat >= 28) & (lat <= 52) & (lon <= -95) & (lon >= -127))
            # 向两边扩一下，这样valid_points给的北美边界的像元能匹配到外部离得最近的点
            indices2 = np.where(valid_points2)  # 获取满足条件的索引

            x_coor2 = x_coor0[indices2]
            y_coor2 = y_coor0[indices2]
            data_select = data[indices2]
            # num2 = int(x_coor2.shape[0])  # 获取矩阵大小，转换为整形。二维  [mm, nn] = np.shape(x_coor)
            # 这里的coordinates1和2对应的都是1500*2500中的索引
            # row_indices2, col_indices2 = indices2

            num1 = 193
            x_input_site = np.zeros((num1, 9), dtype=float)  # 每行对应北美每个测站临近9个像元的x坐标
            y_input_site = np.zeros((num1, 9), dtype=float)  # 每行对应北美每个测站临近9个像元的y坐标
            x_input_grid = np.zeros((382, 9), dtype=float)  # 每行对应北美每个格网点临近9个像元的x坐标
            y_input_grid = np.zeros((382, 9), dtype=float)  # 每行对应北美每个格网点临近9个像元的y坐标

            for iflag1 in range(0, 193):  # 逐个像元遍历
                # print(iflag1)

                # 这里欧几里得距离ss一共有num1个，每个是大小为num2的数组  ss在每一次循环中会更新
                ss_site = np.sqrt(np.square(x_coor2 - x_site[iflag1]) + np.square(y_coor2 - y_site[iflag1]))
                CMI_ind = np.argpartition(ss_site, 9)[:9]  # 提取最小的9个值对应的索引
                ma_ss = 11000
                # 判断 CMI_ind 是否取到,这里是排除整张影像在西北美一个数据都没有的情况（与前面的去掩膜处理对应）
                if len(CMI_ind) == 9:
                    ma_ss = np.max(ss_site[CMI_ind])
                if ma_ss < 10000:  # 确保所有的像元离目标点的距离不超过10km, 加上前面的坐标去掩膜的操作就可以处理那些缺失像元比较多的图像
                    # 这个对应所有点寻找极值
                    x_input_site[iflag1, :] = x_coor2[CMI_ind]  # 北美区每个测站周围临近的9个像元
                    y_input_site[iflag1, :] = y_coor2[CMI_ind]
                    CMI_bt = data_select[CMI_ind]
                    function4 = OrdinaryKriging(x_input_site[iflag1, :], y_input_site[iflag1, :], CMI_bt,
                                                variogram_model='power',
                                                verbose=False, enable_plotting=False,
                                                variogram_parameters=[0.1, 3, 3])
                    CMI_bt_site[irs * 193 + iflag1, 0:2] = (hour2, min2)
                    CMI_bt_site[irs * 193 + iflag1, index_band3 + 2], std_interp = (
                        function4.execute('points', x_site[iflag1], y_site[iflag1]))
                else:
                    stop = 1

            for iflag2 in range(0, 382):  # 逐个像元遍历
                # print(iflag2)

                # 这里欧几里得距离ss一共有num1个，每个是大小为num2的数组  ss在每一次循环中会更新
                ss_site2 = np.sqrt(np.square(x_coor2 - x_grid[iflag2]) + np.square(y_coor2 - y_grid[iflag2]))
                CMI_ind2 = np.argpartition(ss_site2, 9)[:9]  # 提取最小的9个值对应的索引
                ma_ss2 = 11000
                if len(CMI_ind2) == 9:
                    ma_ss2 = np.max(ss_site2[CMI_ind2])
                if ma_ss2 < 10000:  # 确保所有的像元离目标点的距离不超过10km, 加上前面的坐标去掩膜的操作就可以处理那些缺失像元比较多的图像
                    # 这个对应所有点寻找极值
                    # 这里得额外再加一个判断，部分10km范围内的9个最邻近像元缺失，致使结果为9999.9，所以最好是剔除这个历元或者选择临近9个不为空，且不超过10km的历元
                    x_input_grid[iflag2, :] = x_coor2[CMI_ind2]  # 北美区每个测站周围临近的9个像元
                    y_input_grid[iflag2, :] = y_coor2[CMI_ind2]
                    CMI_bt2 = data_select[CMI_ind2]
                    function42 = OrdinaryKriging(x_input_grid[iflag2, :], y_input_grid[iflag2, :], CMI_bt2,
                                                 variogram_model='power',
                                                 verbose=False, enable_plotting=False,
                                                 variogram_parameters=[0.1, 3, 3])
                    CMI_bt_grid[irs * 382 + iflag2, 0:2] = (hour2, min2)
                    CMI_bt_grid[irs * 382 + iflag2, index_band3 + 2], std_interp2 = (
                        function42.execute('points', x_grid[iflag2], y_grid[iflag2]))
                else:
                    stop2 = 1  # 每个历元初始化为0
            # 对于CMI_bt_site、如果是符合10km要求，则正常内插，不然就等于初始化的0 # 在最后一个通道数据赋予完毕的时候判断
            # 对于单历元3通道数据，只要有一个站/格网周围的像元的最大距离超过10 km，则3组数据均赋为0
            if (stop == 1 or stop2 == 1) and index_band3 == 2:
                CMI_bt_site[(irs * 193):(irs * 193 + 193), :] *= 0  # 左闭右开
                CMI_bt_grid[(irs * 382):(irs * 382 + 382), :] *= 0
    # 如果任意
    mask_zero = np.any(CMI_bt_site[:, 2:5] == 0.0, axis=1)  # 判断后三列是否有0
    mask_zero2 = np.any(CMI_bt_grid[:, 2:5] == 0.0, axis=1)
    filtered_CMI_bt_site = CMI_bt_site[~mask_zero]  # 保留非0行，即为所有格网/测站均有内插值的历元
    filtered_CMI_bt_grid = CMI_bt_grid[~mask_zero2]  # 通过最前面+这部分操作，可以确保剔除了缺失通道、以及某一通道数据覆盖范围不全的历元
    np.savetxt("B:\\GOES-16\\result\\CMI_model_" + "{:03}".format(idoy) + ".txt", filtered_CMI_bt_site, fmt='%.5f')
    np.savetxt("B:\\GOES-16\\result\\CMI_grid_" + "{:03}".format(idoy) + ".txt", filtered_CMI_bt_grid, fmt='%.5f')

# Main function to initiate multiprocessing
if __name__ == "__main__":
    pathroot0 = "H:\\GOES-16\\test\\"
    for idoy0 in range(190, 191):
        pathroot2 = pathroot0 + "{:03}".format(idoy0) + "\\"
        process_file(pathroot2, idoy0)
"""
    # Create a Pool of 24 processes
    with Pool(processes=24) as pool:
        pool.map(process_file, era5_single_files)
"""
# 初步可视化

"""
fig0, axes0 = plt.subplots(1, 1, figsize=(4, 3), dpi=300, facecolor="w")
plt.plot(
    blh_num[:, 0],  # X轴数据
    CMI_bt_grid[:, 0],  # Y轴数据
    marker='s',  # 标记样式：圆点
    linestyle='None',  # 线条样式：实线
    color=(169/255, 221/255, 213/255),  # 线条颜色：蓝色
    # linewidth=2,  # 线宽：2
    markersize=6,  # 标记大小：8
    label='Data points (T2m)'  # 图例标签
)

axes0.set_xlabel('Longitude (°W)')
axes0.set_ylabel('CMI_BT (K)')
plt.show()
"""

# 屏蔽无效数据点
"""
data_masked = np.ma.masked_where((data0 == None) | (data0 == '') | np.isnan(data0), data0)
lat_masked = np.ma.masked_where(np.isnan(lat0), lat0)
lon_masked = np.ma.masked_where(np.isnan(lon0), lon0)
"""
