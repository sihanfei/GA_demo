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
    m=m+Pf[i]
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
V = np.zeros(Mn) # 计算值
Pf = np.zeros(Mn) # 归一化生存率
# 产生一个族群
for i in range(Mn):
  G0.append(InitDNA(Bn))
  x = DNA2Data(G0[i]) 
  V[i] = x*np.sin(10.0*np.pi*x)+2.0 # 

Pf = np.divide(V, sum(V))
sortnum = np.argsort(Pf) # 排序对应的序号
Pf = np.sort(Pf) # 生存率排序

# 开始进行遗传变异
Pc = 0.5 # 交叉概率
Pm = 0.05 # 变异概率
T = 500 # 变异代数
E = 0.001 # 停止条件，同代内的最大最小值差

for j in range(T):
  G=[] # 新族群  
  for i in range(Mn/2):
    i0 = RatioSelect(Pf) # 选择DNA序号
    i1 = RatioSelect(Pf)
  
    r = np.random.random() # 交叉
    if r<Pc: 
      DNA0, DNA1 = CrossTrans(G0[sortnum[i0]], G0[sortnum[i1]], Bn)
    else:
      DNA0 = G0[sortnum[i0]]
      DNA1 = G0[sortnum[i1]]
    
    r = np.random.random() # 变异1
    if r<Pm:
      DNA0 = GeneticMutation(DNA0, Bn)
    else:
      DNA0 = DNA0
    
    r = np.random.random() # 变异2
    if r<Pm:
      DNA1 = GeneticMutation(DNA1, Bn)
    else:
      DNA1 = DNA1
    
    if i==0:
      DNA0 = G0[sortnum[-1]] # 保留最好的结果
    
    G.append(DNA0) 
    x = DNA2Data(DNA0) # 计算生存率
    V[i*2] = x*np.sin(10.0*np.pi*x)+2.0
    G.append(DNA1)
    x = DNA2Data(DNA1)
    V[i*2+1] = x*np.sin(10.0*np.pi*x)+2.0
    
  G0 = G # 新族群替代老族群
  Pf = np.divide(V, sum(V))
  sortnum = np.argsort(Pf)
  Pf = np.sort(Pf) # 归一化并排序
  diffV = np.max(V) - np.min(V)
  print '差值为%f' %(diffV)
  if (diffV)<E: 
    if j<150:
      Pm = 0.1 # 提高变异率
    else:
      break # 达到结束条件
  
x = DNA2Data(G0[sortnum[-1]]) 
print '遗传代数为%d, 最大值为%f,最小值为%f,此时的x为%f' %(j, np.max(V), np.min(V), x)
