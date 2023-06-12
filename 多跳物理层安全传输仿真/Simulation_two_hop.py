import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd


def ppp_2d(am, sign, lam):
    r = 10       # 假设平面的范围为方圆10km
    pi = math.pi
    lam2 = am * lam * pi * pow(r, 2)     # 新的密度（结合面积而来）
    n = np.random.poisson(lam=lam2, size=1)  # 产生泊松点的数目
    u1 = np.sort(np.random.uniform(0, 1, n))  # 均匀分布
    R = [r * math.sqrt(i) for i in u1]  # 距离
    u2 = np.random.uniform(0, 1, n)
    W = [2 * pi * j for j in u2]  # 角度

    if sign == 1:                 # 目的是辨别基站是否含有请求内容
        Sign = np.ones(n)
    else:
        Sign = np.zeros(n)

    Loc = list(zip(R, W, Sign))

    #print(Sign)
    #print(u1)
    #print(R)
    #print(u2)
    #print(W)
    #print(Loc)
    return Loc


def distance(a, b):      # 两个点之间的距离计算，输入要求是一个（和原点的距离，角度）
    return math.sqrt(pow(a[0], 2)+pow(b[0], 2)-2*a[0]*b[0]*math.cos(a[1]-b[1]))


def max_distance(p=1000, Rs=5, Rv=18, sigma=1, a=4):  # 基站的覆盖范围
    return pow(p/(pow(2, Rs+Rv)*sigma), 1/a)


def metric(p=1000, Rv=18, sigma=1):     # 窃听成功：power(r1,-a)+power(r2,-a) >= metric
    return (pow(2, Rv)*sigma)/p


def first_hop(loc1, loc2):                   # 判断合法用户是否可以成功连接到基站
    d = max_distance()
    if loc2:
        d2 = loc1[0][0]                         # 由于我们之前排了个序，所以loc[0]肯定是距离原点最近的点
        d3 = loc2[0][0]
        loc = loc1[0] if d2 < d3 else loc2[0]   # 选择最近的基站
    else:
        loc=loc1[0]
    if loc[0] < d:
        if loc[2] == 1:
            return 1, loc                      # 成功连接，且基站有所请求内容,loc是基站信息
        else:
            return 2, loc                       # 成功连接，但基站没有所请求内容,loc是基站信息
    else:
        return 0,loc                          # 连接失败


def eves_first_hop(loc, loc3, a=4):   # 针对连接基站的联合窃听
    metric_v = metric()                # 确认指标，这个指标只与Rv相关
    d = []
    for i in loc3:
        d1 = distance(i, loc)
        d.append(d1)
    d.sort()                           # 计算所有窃听者与第一跳基站的距离，并排序
    # print(d)

    if pow(d[0], -a) + pow(d[1], -a) <= metric_v:  # 判断联合窃听成功与否
        # if pow(d[0], -a) <= metric_v:
        return 0                                   # 窃听失败
    else:
        return 1                                   # 窃听成功

'''
  只有在第一跳成功连接，但是没有所请求的内容，同时第一跳的联合窃听失败的情况下，才会尝试连接第二跳。
'''


def second_hop(loc1, loc):                  # 第二跳是有所请求的内容且距离第一跳最近的基站
    d = max_distance()
    dis = []
    for i in loc1:
        d1 = distance(i, loc)
        dis.append((d1,i))
    dis = sorted(dis,key=lambda x:x[0])     # 计算所有有请求内容的基站与第一跳基站的距离，并排序
    # print(dis)
    loc2 = dis[0][1]                        # 距离最近的基站位置
    if dis[0][0] < d:
        return 1, loc2                      # 成功连接，且基站有所请求内容,loc是基站信息
    else:
        return 0, loc2                           # 连接失败


def eves_second_hop(loc2, loc3, a=4):         # 针对连接基站的联合窃听
    metric_v = metric()                       # 确认指标，这个指标只与Rv相关
    d = []
    for i in loc3:
        d1 = distance(i, loc2)
        d.append(d1)
    d.sort()                                  # 计算所有窃听者与第一跳基站的距离，并排序
    # print(d)

    if pow(d[0], -a) + pow(d[1], -a) <= metric_v:    # 判断联合窃听成功与否
        # if pow(d[0], -a) <= metric_v:
        return 1                                     # 窃听失败
    else:
        return 0                                     # 窃听成功


f_res = 0


def start(Am):                                     # Am是文件存储概率
    global f_res
    data = []
    s = np.arange(10, 100, 10)
    print(s)
    for lamb in s:
        sum1 = 0
        i = 0
        j = 0
        Am = Am
        m = np.random.zipf(a=2, size=1000)         # 产生1000个数字的满足zipf分布的序列
        print(m)
        while i < 1000:
            lame = 1
            if m[i] < len(Am):
                j = j + 1
                am = Am[m[i]-1]

                if am == 0:
                    f_res = 0
                elif am == 1:
                    loc1_all = ppp_2d(am=am, sign=1, lam=lamb)            # 含有请求内容的基站
                    loc3_all = ppp_2d(am=1, sign=0, lam=lame)           # 窃听者
                    res1, loc1 = first_hop(loc1_all, None)
                    if res1 == 0:                               # 如果第一跳未能建立连接
                        f_res = 0                               # 判定：失败
                    elif res1 == 1:                             # 如果第一跳建立了连接且有请求内容
                        res2 = eves_first_hop(loc1, loc3_all, a=4)  # 判断是否窃听成功
                        if res2 == 0:                           # 窃听失败
                            f_res = 1                           # 判定：成功
                        elif res2 == 1:                         # 窃听成功
                            f_res = 0                           # 判定：失败
                else:
                    loc1_all = ppp_2d(am=am, sign=1, lam=lamb)         # 含有请求内容的基站
                    loc2_all = ppp_2d(am=1-am, sign=0, lam=lamb)       # 不含请求内容的基站
                    loc3_all = ppp_2d(am=1, sign=0, lam=lame)        # 窃听者
                    res1,loc1 = first_hop(loc1_all,loc2_all)
                    if res1 == 0:                             # 如果第一跳未能建立连接
                        f_res = 0                             # 判定：失败
                    elif res1 == 1:                           # 如果第一跳建立了连接且有请求内容
                        res2 = eves_first_hop(loc1, loc3_all, a=4) # 判断是否窃听成功
                        if res2 == 0:                           # 窃听失败
                            f_res = 1                         # 判定：成功
                        elif res2 == 1:                         # 窃听成功
                            f_res = 0                         # 判定：失败
                    elif res1 == 2:                           # 如果第一跳建立了连接且没有请求内容
                        res2 = eves_first_hop(loc1, loc3_all, a=4) # 判断是否窃听成功
                        if res2 == 1:                           # 窃听成功
                            f_res = 0                         # 判定：失败
                        elif res2 == 0:                         # 窃听失败
                            res3, loc2 = second_hop(loc1_all, loc1)      # 尝试连接第二跳基站
                            if res3 == 0:                       # 未能找到可以含有请求建立连接的第二跳
                                f_res = 0                     # 判定：失败
                            elif res3 == 1:                     # 找到可以含有请求建立连接的第二跳
                                f_res = eves_second_hop(loc2, loc3_all, a=4) # 判定是否窃听成功
                sum1 = sum1 + f_res
            i = i + 1
        data.append(sum1/j)
        print(lamb)
        print(data)
    return s, data

# start(Am=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1])


s, data = start(Am=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.1])
data2 = {'index': s, 'data': data}
pd.DataFrame(data2).to_csv(r'G:\\data_lam_1.2.csv')












