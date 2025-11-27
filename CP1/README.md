Codes for **Computer Project 1** of course 'Applied Stochastic Analysis' in Peking University. For details see:

https://www.math.pku.edu.cn/teachers/litj/notes/appl_stoch/ComputerProjects2024.pdf

# Discription of codes

[Potts_Model_2D.py](Potts_Model_2D.py): 实现了PottsModel2D类,  这个类包含了Potts Model的各个参数, 以及储存了当前状态下的lattice. 这个类的方法包括: 初始化各个参数, 设置温度, 计算单个状态下的物理量, 实现一次Wolff算法, 实现一次single-flip算法等. 此外, 这个文件还实现了MCMC模拟的功能, 即给定Potts Model的各个参数和各阶段迭代次数, 返回需要的物理量.