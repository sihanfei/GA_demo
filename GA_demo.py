#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 10:43:25 2017

@author: taosheng
"""

import numpy as np
from bitarray import bitarray

"""############################################################################
purpse：随机选择1个DNA
input 
  Pf: 各DNA生存率（按照升序排列）
output：
  i：被选中的DNA的序号
"""
def RatioSelect(Pf):
  m=0
  r=np.random.random()
  for i in range(np.size(Pf)):
    m=m+r
    if r<m:
      return i


"""############################################################################
purpse：选择一个随机位置进行两个DNA的粒交换
input： 
  M1，M2：输入DNA
  Bn 总的位数
output：
  Cm1，Cm2：交换后的DNA
"""
def CrossTrans(M1, M2, Bn):
  Cm1 = bitarray(Bn)
  Cm2 = bitarray(Bn)
  p = np.random.randint(Bn-1) # 选择交换位
  Cm1[:p]=M1[:p] # 交换
  Cm1[p:]=M2[p:]

  Cm2[:p]=M2[:p]
  Cm2[p:]=M1[p:]   
  return Cm1, Cm2    

"""############################################################################
purpse：选择一个随机位置进行DNA的变异
input： 
  M：输入DNA
  Bn：总的位数
output：
  M：变换后的DNA
"""
def GeneticMutation(M, Bn):
  p = np.random.randint(Bn-1) # 选择变异位
  M[p] ^= 1 # 变异
  return M

"""############################################################################
purpose：根据DNA的粒子数，生成一个初始DNA
input： 
  Bn：DNA的粒子数，也就是位数，和问题直接相关
output：
  Mi：初始化的DNA
"""
def InitDNA(Bn):
  Ri = np.random.randint(2**Bn) # 产生一个随机数 0～2^Bn
  SRi = str(bin(Ri)) # 转换为str
  BAi = bitarray(SRi[2:]) # 转换为bitarray    
  Mi = bitarray(Bn) # 生成空的bitarray
  Mi[0:BAi.length()]=BAi
  return Mi

"""############################################################################
purpose：从bitarray转换成int
input： 
  BA：单个bitarray
output：
  INT：转换完的结果
"""
def fromBA2INT(BA):
  S = BA.to01()
  INT = int('0b'+S, 2)
  return INT

"""############################################################################
遗传算法实际：f(x)=x*sin(10*pi*x)+2, x~[-1,2],精度0.01
Bn：（2-（-1））／0.01=300<2^9
"""

"""############################################################################
purpose：从DNA转到实际数据
input： 
  DNA：单个bitarray
output：
  Data：实际数据（float）
"""
def DNA2Data(DNA):
  value = fromBA2INT(DNA)
  return (value*3.0/512.0-1)

Mn=20 # 20个DNA
Bn=9 # DNA有9个粒
G0=[] # 族群
Pf=np.zeros(Mn) # 生存率
# 产生一个族群
for i in range(Mn):
  G0.append(InitDNA(Bn))
  x = DNA2Data(G0[i])
  Pf[i] = x*np.sin(10.0*np.pi*x)+2.0

Pf = np.divide(Pf, sum(Pf))
Pf= np.sort(Pf) # 归一化并排序

# 
