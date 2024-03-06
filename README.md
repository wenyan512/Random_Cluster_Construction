# Random Cluster Construction Model

## 概述

Packmol是一款构建复杂体系的软件，通过配置文件，可在指定空间区域中堆积分子，创建分子系统模型。详见：https://m3g.github.io/packmol/。 本工作旨在开发一个Python程序（脚本），构建随机团簇（cluster）模型。

Packmol的用法可见：

1. https://m3g.github.io/packmol/examples.shtml

2. https://blog.chembiosim.com/Packmol-basic-01/

## 功能

本程序主要功能为：给定组分分子结构文件以及指定建模参数，调用Packmol程序生成最优的随机团簇模型，并输出相关信息。

具体要求如下：

1. 支持球形或方形模型；

2. 支持单组分与多组分；

3. 指定各组分分子数量；

4. 指定尺寸或自动计算最佳尺寸；

5. 返回模型文件和盒子尺寸。

## 用法

输入文件：一个或多个组分分子文件（mol2格式）

建模参数：

| 名称         | 类型        | 含义 
| ----------- | ----------- | ----------- |
| shape       | enum[“sphere”, “cuboid”]            | 模型形状 |
| components  | [{“file”: str, “count”: int}, ...]  | 组分分子文件及对应数量 |
| size        | (x, y, z) \| None                    |盒子尺寸。若为None，则自动计算最佳尺寸。|

输出数据：

| 名称       | 类型               | 含义 
| --------- | -----------        | ----------- |
| file      | str                | 模型文件（pdb格式，可固定为model.pdb）|
| size      | (x, y, z) \| None   | 盒子尺寸。若为None，则自动计算最佳尺寸。|
