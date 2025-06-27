# TBkit

parkman@buaa.edu.cn

fanyang_tcmp@buaa.edu.cn

yzpeng@buaa.edu.cn

You_YiQun@buaa.edu.cn


#### 介绍
TBkit是 BUAA-TCMP 推出的基于 matlab 的第一性原理计算结果分析处理，紧束缚模型与低能有效kp模型的构建，调控与计算。为拓扑材料，磁性材料以及量子输运相关科研工作者提供稳定快速、简易高效的程序开发框架
无论是个人、团队、或是高校，都能够用 TBkit实...。

如有什么问题可以挂issue或直接发邮箱；欢迎您的建议

#### 软件架构
外部程序：spglib
matlab工具箱：
TBkit基本函数，TBkit类，HR类，HK类，Htrig类


![alt text](doc/TBKit_.png)


#### 安装教程
运行INSTALL.m脚本

run the INSTALL.m script

#### 使用说明

直接使用；编写脚本或命令行使用

#### 短期目标
- [ ] Clifford 类的建立 一般性地对哈密顿量进行Clifford拆分
- [ ] 强化 HR Htrig HK 对称性分析功能
- [ ] 向量化 对称性 Htrig HK
- [ ] 异质结
- [ ] Hckt针对一般TB(带复数)的转化 与实验数据的对接
- [ ] IBZ

#### 发展目标

- [ ] TB 模型中各种场的加入，拉伸以及缺陷，散射！
- [ ] Oper Class 的特征标构建 （waiting）
- [ ] HCEF的构建 
- [ ] 构建 各类的 机器学习拟合DFT能带的程序 （2021.07.08-07.12）完成80%
- [ ] （待选）学习牛派的半经典方法，后续接入VASP的WAVECAR读取，先做个平面波的方法，再做个VASP的借口
- [ ] （待选）神经网络的matlab实现，对接DFT或拓扑以及拓扑量子化学 （）
- [ ] （待选）自旋模型的对称性构建，自动生成 简单ED模型求解基本强关联特性



