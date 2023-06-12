import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd


def ppp_2d(am, lam):
    r = 10       # 假设平面的范围为方圆10km
    pi = math.pi
    lam2 = am * lam * pi * pow(r, 2)     # 新的密度（结合面积而来）
    n = np.random.poisson(lam=lam2, size=1)  # 产生泊松点的数目
    u1 = np.sort(np.random.uniform(0, 1, n))  # 均匀分布
    R = [r * math.sqrt(i) for i in u1]  # 距离
    u2 = np.random.uniform(0, 1, n)
    W = [2 * pi * j for j in u2]  # 角度
    '''
    if sign == 1:
        Sign = np.ones(n)
    else:
        Sign = np.zeros(n)
    '''
    Loc = list(zip(R, W))

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


def first_hop(loc1):                   # 判断合法用户是否可以成功连接到基站
    d = max_distance()
    d2 = loc1[0][0]                    # 由于我们之前排了个序，所以loc[0]肯定是距离原点最近的点

    if d2 < d:
        return 1                       # 成功
    else:
        return 2                       # 失败


def eves_first_hop(loc1, loc3, a=4):   # 针对连接基站的联合窃听
    metric_v = metric()                # 确认指标，这个指标只与Rv相关
    d = []
    for i in loc3:
        d1 = distance(i, loc1[0])
        d.append(d1)
    d.sort()                           # 计算所有窃听者与基站的距离，并排序
    # print(d)

    if pow(d[0], -a) + pow(d[1], -a) <= metric_v:  # 判断联合窃听成功与否
    #if pow(d[0], -a) <= metric_v:
        return 1                                   # 成功
    else:
        return 0                                   # 失败


def start(Am):                                     # Am是文件存储概率
    data = []
    s = np.arange(1, 11, 1)
    print(s)
    for lame in s:
        sum1 = 0
        i = 0
        j = 0
        Am = Am
        m = np.random.zipf(a=2, size=10000)         # 产生1000个数字的满足zipf分布的序列
        print(m)
        while i < 10000:

            if m[i] < len(Am):
                j = j + 1
                am = Am[m[i]-1]

                if am == 0:
                    f_res = 0

                else:
                    loc1 = ppp_2d(am=am, lam=55)
                    loc3 = ppp_2d(am=1, lam=lame)
                    res1 = first_hop(loc1)
                    if res1 == 1:
                        res2 = eves_first_hop(loc1, loc3, a=4)
                        #res2 = 1
                    elif res1 == 2:
                        res2 = 0
                    f_res = res2

                sum1 = sum1 + f_res
            i = i + 1
        data.append(sum1/j)
        print(lame)
        print(data)

#start(Am=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
start(Am=[1.0000, 1.0000, 0.7373, 0.5691, 0.4573, 0.3711, 0.3002, 0.2395, 0.1864, 0.1390])













