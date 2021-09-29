## 微分求积法（Differential Quadrature）学习笔记

## 1.DQ方法的基本原理

微分求积法(DQ)是Bellman和他的同事在70年代初提出的。DQ方法是求导数近似的一种数值离散化方法，起源于传统积分求积的思想。

### 1.1 积分求积

一般地，积分$\int_a^bf(x)dx$可以近似为：
$$
\int_a^bf(x)dx = w_1f_1+w_2f_2+\dots + w_nf_n = \sum_{k=1}^n w_k f_k
\tag{1.1}
\label{1.1}
$$
其中$w_1,w_2,\dots,w_n$为权值系数，$f_1,f_2,\dots,f_n$为函数在离散点$a=x_1,x_2,\dots,x_n = b$处的值。方程$(\ref{1.1})$称为积分求积，它使用整个积分域中的所有函数值来近似一个有限积分上的积分。

### 1.2 微分求积

函数$ f (x)$对 $x $在网格点$x_i$处的一阶导数可以近似为整个域中所有函数值的线性和：
$$
f_x(x_i)=\frac{df}{dx}\vert_{x_i} = \sum_{j=1}^N a_{ij}\cdot f(x_j),~~~i=1,2,\dots,N
\tag{1.2}
\label{1.2}
$$
其中$ a_{ij}$ 表示加权系数，$N$表示整个域中的网格点数。方程$(\ref{1.2})$称为微分求积(DQ).

## 2.基于多项式的微分求积法(PDQ)

在DQ近似中，如何确定加权系数很多工作都是基于多项式逼近的，因此相关的方法可以被称为基于多项式的微分求积法(PDQ)。

### 2.1  一阶导数加权系数的计算

Bellman 等人假定函数$ f (x)$在区间$[ a,b ]$上足够光滑，因此它在任何网格点上的一阶导数 $f^{ (1)}(x)$可以用下面的公式来近似：
$$
f^{(1)}_x(x_i)=\frac{df}{dx}\vert_{x_i} = \sum_{j=1}^N a_{ij}\cdot f(x_j),~~~i=1,2,\dots,N
\tag{2.1}
\label{2.1}
$$
方程$(\ref{2.1})$中加权系数$a_{ij}$的确定是DQ近似的关键步骤。下面是一些方法

#### 2.1.1 Bellman方法

- **第一个方法**.   使用测试函数如下：

$$
g_k(x)=x^k,~~k=0,1,\dots,N-1
\tag{2.2}
\label{2.2}
$$

显然有$N$个测试函数，那么需要N个网格点$x_1,x_2,\dots,x_N$,可以得到关于$a_{ij}$的$N \times N$个线性代数方程如下：
$$
\left\{\begin{array}{l}
\sum_{j=1}^{N} a_{i j}=0 \\
\sum_{j=1}^{N} a_{i j} \cdot x_{j}=1 \\
\sum_{j=1}^{N} a_{i j} \cdot x_{j}^{k}=k \cdot x_{i}^{k-1}, k=2,3, \ldots, N-1
\end{array}\right.~~~~i=1,2, \cdots, N.
\tag{2.3}
\label{2.3}
$$
<span style='color:greenyellow'>由$\color{greenyellow} k=0,k=1,k=2,3, \ldots, N-1$分别代入得到。</span>

方程组2.3有唯一解，因为它的矩阵是范德蒙形式：
$$
{\rm{V = }}\left[ {\begin{array}{*{20}{c}}
{\rm{1}}&{\rm{1}}&{\rm{1}}& \cdots &{\rm{1}}\\
{{x_1}}&{{x_2}}&{{x_3}}& \cdots &{{x_N}}\\
{x_1^2}&{x_2^2}&{x_3^2}& \cdots &{x_N^2}\\
 \vdots & \vdots & \vdots & \ddots & \vdots \\
{x_1^{N - 1}}&{x_2^{N - 1}}&{x_3^{N - 1}}& \cdots &{x_N^{N - 1}}
\end{array}} \right]
\notag
$$


但是，当$n$很大的时候，矩阵是病态的，求解它的逆是困难的。在这种方法的实际应用中，通常选择$\color{orchid}n<13$。

- **第二个方法**.  与第一种方法类似，但是使用不同的测试函数，为：
  $$
  g_k(x) = \frac{L_N(x)}{(x-x_k)\cdot L_N^{(1)}(x_k)},~~k = 1,2,\dots,N
  \tag{2.4}
  \label{2.4}
  $$
  其中$L_N(x)$为$N$次Legendre多项式，$L_N^{(1)}$为其一阶导数。这里需选择$x_k$为shifted Legendre多项式的根.在$N$个网格点上应用$(\ref{2.4})$得到关于$a_{ij}$的代数公式
  $$
  a_{i j}=\frac{L_{N}^{(1)}\left(x_{i}\right)}{\left(x_{i}-x_{j}\right) L_{N}^{(1)}\left(x_{j}\right)}, \quad \text { for } j \neq i 
  \tag{2.5a}
  \label{2.5a}
  $$

  $$
  \color{brown}
  a_{i i}=\frac{1-2 x_{i}}{2 x_{i}\left(x_{i}-1\right)}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  \tag{2.5b}
  \label{2.5b}
  $$

  <span style='color:greenyellow'>由$\color{greenyellow}{(\ref{2.4})}$可以得到，$\color{greenyellow}g_k(x_j)=\delta_{kj}$,那么$\color{greenyellow}k \neq i$时：</span>

$$
\color{greenyellow}
g_k'(x_i) = \sum_{j=1}^N a_{ij} g_k(x_j)=a_{ij},k=j，\\
\color{greenyellow}
g_k’(x) = \frac{L_N^{(1)}(x)(x-x_k)-L_N(x)}{(x-x_k)^2L_N^{(1)}(x_k)},\\
\color{greenyellow}
g_k'(x_i) =\frac{L_{N}^{(1)}\left(x_{i}\right)}{\left(x_{i}-x_{j}\right) L_{N}^{(1)}\left(x_{j}\right)}
\notag
$$

<span style='color:greenyellow'> $\color{greenyellow}k = i$时，根据L‘Hospital规则以及$\color{greenyellow}P_N(x)$满足的条件</span>
$$
\color{greenyellow}
(1-x^2)L_N''(x)-2xL_N'(x)+N(N+1)L_N(x)=0
\notag
$$
<span style='color:greenyellow'>得到：</span>
$$
\color{greenyellow}
a_{ii}=\lim_{x \to xi}\frac{L_N^{(1)}(x)(x-x_i)-L_N(x)}{(x-x_i)^2L_N^{(1)}(x_i)}\\
\color{greenyellow}
=\lim_{x \to xi}\frac{L_N^{(2)}(x)(x-x_i)+L_N^{(1)}(x)-L_N^{(1)}(x)}{2(x-x_i)L_N^{(1)}(x_i)}\\
\color{greenyellow}
=\lim_{x \to xi}\frac{L_N^{(2)}(x)}{2L_N^{(1)}(x_i)}=\lim_{x \to xi}\frac{2xL_N'(x)-N(N+1)L_N(x)}{(1-x^2)L_N^{(1)}(x_i)}\\
\color{greenyellow}
=\frac{2x_i}{2(1-x_i^2)}
\notag
$$
这种方法不像第一种方法那样灵活，因为这种方法中网格点的坐标不能任意选择，所以在实际应用中通常采用第一种方法。

#### 2.1.2 Quan and Chang's 方法

Quan and Chang用下面的拉格朗日插值多项式作为检验函数：
$$
g_k(x)=\frac{M(x)}{(x-x_k)\cdot M^{(1)}(x_k)},~~k=1,2,\dots,N
\tag{2.6}
\label{2.6}
$$
其中:
$$
M(x)=(x-x_1)(x-x_2)\dots(x-x_N)
\tag{2.7}
\label{2.7}
$$

$$
M^{(1)}(x_i)= \prod \limits_{k=1,k\neq i}^{N}(x_i-x_k)
\tag{2.8}
\label{2.8}
$$

随后，通过在$N$个网格点上应用方程$(\ref{2.6})$，得到了以下计算加权系数$a_{ij}$的代数公式：
$$
a_{i j}=\frac{1}{x_{j}-x_{i}} \prod_{k=1, k \neq i, j}^{N} \frac{x_{i}-x_{k}}{x_{j}-x_{k}}, \text { for } j \neq i
\tag{2.9a}
\label{2.9a}
$$

$$
a_{i i}=\sum_{k=1, k \neq i}^{N} \frac{1}{x_{i}-x_{k}}
\tag{2.9b}
\label{2.9b}
$$

显然可以得到$g_k(x_j)=\delta_{jk}$
$$
\color{greenyellow}
g_k'(x_i) = \sum_{j=1}^N a_{ij}g_k(x_j)=a_{ik}\\
\color{greenyellow}
g_k'(x) = \frac{M'(x)(x-x_k)-M(x)}{(x-x_k)^2\cdot M^{(1)}(x_k)}\\


\notag
$$
<span style='color:greenyellow'>当$\color{greenyellow} i \neq j $时</span>
$$
\color{greenyellow}
\therefore~~ g_k'(x_i) = \frac{M'(x_i)}{(x_i-x_k)\cdot M^{(1)}(x_k)}=\frac{M^{(1)}(x_i)}{(x_i-x_k)\cdot M^{(1)}(x_k)}\\
\color{greenyellow}
\therefore ~~a_{ik}=\frac{M^{(1)}(x_i)}{(x_i-x_k)\cdot M^{(1)}(x_k)} = \frac{1}{(x_i-x_k)}\cdot \frac{\prod \limits_{m=1,m\neq i}^N (x_i-x_m)}{\prod \limits_{m=1,m\neq k}^N (x_k-x_m)}\\
\color{greenyellow}
\therefore ~~a_{ij}=\frac{1}{x_{j}-x_{i}} \prod_{k=1, k \neq i, j}^{N} \frac{x_{i}-x_{k}}{x_{j}-x_{k}}
\notag
$$
<span style='color:greenyellow'>当$\color{greenyellow} i =j$ 时</span>
$$
\color{greenyellow}
a_{ii} = \lim_{x \to x_i}\frac{M'(x)(x-x_i)-M(x)}{(x-x_i)^2M^{(1)}(x_i)}\\
\color{greenyellow}
= \lim_{x \to x_i}\frac{M''(x)}{2M^{(1)}(x_i)}=\lim_{x \to x_i}\frac{2\sum\limits_{m=1,m\neq i}^N\prod\limits_{k=1,k\neq i, m}^N(x-x_k)}{2\prod_{k=1,k\neq i}^N(x_i-x_k)}\\
\color{greenyellow}
=\sum_{k=1,k\neq i}^N \frac{1}{x_i-x_k}
\notag
$$

另一种方法(也就是下面的$(\ref{2.10c})$)：令$M(x)=N(x,x_k)\cdot (x-x_k),k=1,2,\dots,N$，这里$N(x_i,x_j)=M^{(1)}(x_i)\cdot \delta_{ij}$.得到：
$$
g_k(x)=\frac{N(x,x_k)}{M^{(1)}(x_k)},\quad k=1,2,\dots,N\\
\notag
\therefore \quad  a_{ij}=\frac{N^{(1)}(x_i,x_j)}{M^{(1)}(x_j)} \\
$$
又有
$$
M^{(m)}(x)=N{(m)}(x,x_k)\cdot (x-x_k)+m \cdot N{(m-1)}(x,x_k),m=1,2,\dots,N-1;k=1,2,\dots,N 
\tag{2.10}
\label{2.10}
$$

$$
\therefore N^{(1)}(x_i,x_j) = \frac{M^{(1)}(x_i)}{x_i-x_j},i \neq j\\
\quad N^{(1)}(x_i,x_i) = \frac{M^{(2)}(x_i)}{2}\\
\therefore a_{ij} = \frac{M^{(1)}(x_i)}{(x_i-x_j)M^{(1)}(x_j)},i\neq j \\
\quad a_{ii} = \frac{M^{(2)}(x_i)}{2M^{(1)}(x_i)}
\notag
$$

#### 2.1.3  Shu's Approach

偏微分方程的解可以用高次多项式精确逼近。现在。我们假设逼近多项式的次数是$N-1$，这个近似多项式通过向量加法和标量乘法运算构成了一个$N$维线性向量空间，可以用不同的形式表达。基多项式的四个典型集如下:
$$
r_{k}(x)=x^{k-1}, k=1,2, \ldots, N
\tag{2.11a}
\label{2.11a}
$$

$$
r_{k}(x)=\frac{L_{N}(x)}{\left(x-x_{k}\right) \cdot L_{N}^{(1)}\left(x_{k}\right)}, k=1,2, \cdots, N
\tag{2.11b}
\label{2.11b}
$$

$$
r_{k}(x)=\frac{M(x)}{\left(x-x_{k}\right) \cdot M^{(1)}\left(x_{k}\right)}, k=1,2, \cdots, N
\tag{2.11c}
\label{2.11c}
$$

$$
r_{1}(x)=1, r_{k}(x)=\left(x-x_{k-1}\right) \cdot r_{k-1}(x), k=2,3, \ldots, N
\tag{2.11d}
\label{2.11d}
$$

其中$L_N(x)$是N次的Legendre多项式，$m (x)$在方程$(\ref{2.7})$中定义。在四组基多项式中。方程$(\ref{2.11b})$ 和$(\ref{2.11c})$ 来自拉格朗日插值多项式，方程$(\ref{2.11d})$ 来自牛顿插值多项式。方程$(\ref{2.11b})$ 和$(\ref{2.11c})$ 的区别在于网格点的分布。方程$(\ref{2.11b})$是方程$(\ref{2.11c})$的一个特例，因为它只在Legendre配点处有效。

此外，根据线性向量空间的性质，如果一组基多项式满足一个线性算子，那么另一组基多项式也满足。因此 $a_{ij}$满足以下方程，这个方程是由基多项式$x_{k-1}$在$k= 1$时得到的：
$$
\sum_{j=1}^N a_{ij}=0 或 a_{ii} =  - \sum_{j=1,j\neq i}^N a_{ij}
\tag{2.12}
\label{2.12}
$$

### 2.2  二阶导数加权系数的计算

对于二阶导数的离散化，引入了一个类似的近似形式：
$$
f^{(2)}_x(x_i) = \sum_{j=1}^N b_{ij}\cdot f(x_j),~~~i=1,2,\dots,N
\tag{2.13}
\label{2.13}
$$

#### 2.2.1 Quan and Chang's 方法

$$
\left.b_{i j}=\frac{2}{x_{j}-x_{i}}\left(\prod_{k=1, k \neq i, j}^{N} \frac{x_{i}-x_{k}}{x_{j}-x_{k}}\right) \left(\sum_{l=1, l \neq i, j}^{N} \frac{1}{x_{i}-x_{l}}\right)\right., for\quad i \neq j
\tag{2.14a}
\label{2.14a}
$$

$$
b_{i i}=2 \sum_{k=1, k \neq i}^{N-1}\left[\frac{1}{x_{i}-x_{k}}\left(\sum_{l=k+1, l \neq i}^{N} \frac{1}{x_{i}-x_{l}}\right)\right]
\tag{2.14b}
\label{2.14b}
$$

$$
注意：
\color{orange}
g''_k(x_i) = \sum_{j =1}^N b_{ij}g_k(x_j)  = \sum_{m=1}^N a_{im}\sum_{j=1}^N a_{mj} g_k(x_j) \\
\color{orange}
\therefore b_{ij} = \sum_{m=1}^N a_{im}a_{mj} \\
\notag
%\color{greenyellow} 
%\therefore b_{ij}=

%g''(x)=\frac{M''(x)(x-x_k)^2-2M'(x)(x-x_k)+M(x)}{(x-x_k)^3M^{(1)}(x_k)}\\
%\therefore g''_k(x_i) = \frac{M''(x_i)(x_i-x_k)-2M^{(1)}(x_i)}{(x_i-x_k)^2M^{(1)}(x_k)} \\
%\therefore i \neq j时，a_{ij}=g''_k(x_i) = \frac{2\sum \limits_{m=1,m\neq i}^{N} \prod \limits_{h=1,h\neq i,m}^{N}(x_i-x_h)}{(x_i-x_k)M^{(1)}(x_k)}-\frac{2M^{(1)}(x_i)}{(x_i-x_k)^2M^{(1)}(x_k)}\\
%= \frac{2}{x_i-x_k}\left(  \right)
$$

#### 2.2.2 Shu‘s 方法

$$
\because g_k(x)=\frac{N(x,x_k)}{M^{(1)}(x_k)} \\
\therefore b_{ij} = \frac{N^{(2)}(x_i,x_j)}{M^{(1)}(x_j)} \\
又有 N^{(2)}(x_i,x_j) = \frac{M^{(2)}(x_i)-2N^{(1)}(x_i,x_j)}{x_i-x_j} ,i \neq j\\
N^{(2)}(x_i,x_i) = \frac{M^{(3)}(x_i)}{3}\\
\notag
$$

$$
b_{ij} = \frac{M^{(2)}(x_i)-2N^{(1)}(x_i,x_j)}{(x_i-x_j)M^{(1)}(x_j)} ,\quad i\neq j,\\
\tag{2.15a}
\label{2.15a}
$$

$$
b_{ii}= \frac{M^{(3)}(x_j)}{3M^{(1)}(x_j)}
\tag{2.15b}
\label{2.15-b}
$$

由$a_{ij}，a_{ii}$的定义可得：
$$
b_{ij}=2a_{ij}(a_{ii}-\frac{1}{x_i-x_j}),i \neq j
\tag{2.16a}
\label{2.16a}
$$
同样，也可以根据基函数为$x^{k-1}$时的定义得到：
$$
\sum_{j=1}^N b_{ij}=0 或 b_{ii} = -\sum_{j=1,j\neq i}^N b_{ij}
\tag{2.16b}
\label{2.16b}
$$

### 2.3 Shu's方法关于高阶导数的递推公式

$$
f_{x}^{(m-1)}\left(x_{i}\right)=\sum_{j=1}^{N} w_{i j}^{(m-1)} \cdot f\left(x_{j}\right),
\tag{2.17}
\label{2.17}
$$

$$
f_{x}^{(m)}\left(x_{i}\right)=\sum_{j=1}^{N} w_{i j}^{(m)} \cdot f\left(x_{j}\right)
\tag{2.18}
\label{2.18}
$$

其中$i=1,2, \ldots, N ; m=2,3, \ldots, N-1$。同样，我们可以得到：
$$
w_{i j}^{(m-1)}=\frac{N^{(m-1)}\left(x_{i}, x_{j}\right)}{M^{(1)}\left(x_{j}\right)}
\tag{2.19}
\label{2.19}
$$
$$
w_{i j}^{(m)}=\frac{N^{(m)}\left(x_{i}, x_{j}\right)}{M^{(1)}\left(x_{j}\right)}
\tag{2.20}
\label{2.20}
$$
所以得到：
$$
N^{(m-1)}\left(x_{i}, x_{j}\right)=w_{i j}^{(m-1)} M^{(1)}\left(x_{j}\right),i\neq j
\tag{2.21}
\label{2-21}
$$
通过$(\ref{2.10})$，得到:
$$
N^{(m-1)}\left(x_{i}, x_{i}\right)=\frac{M^{(m)}\left(x_{i}\right)}{m}
\tag{2.22}
\label{2.22}
$$

$$
N^{(m)}\left(x_{i}, x_{j}\right)=\frac{M^{(m)}\left(x_{i}\right)-m N^{(m-1)}\left(x_{i}, x_{j}\right)}{x_{i}-x_{j}}, i \neq j \\
N^{(m)}\left(x_{i}, x_{i}\right)=\frac{M^{(m+1)}\left(x_{i}\right)}{m+1}
\tag{2.23}
\label{2.23}
$$

所以：
$$
N^{(m)}(x_i,x_j)=\frac{m[N^{(m)}\left(x_{i}, x_{j}\right)-N^{(m-1)}\left(x_{i}, x_{j}\right)]}{x_i-x_j},i\neq j
\tag{2.24}
\label{2.24}
$$

$$
\therefore~~~N^{(m)}(x_i,x_j)=\frac{m[w_{ii}^{(m-1)}M^{(1)}(x_i)-w_{ij}^{(m-1)}M^{(1)}(x_j)}{x_i,x_j},i \neq j
\tag{2.25}
\label{2.25}
$$

从而有：
$$
w_{ij}^{(m)}=m\left( a_{ij}w_{ii}^{(m-1)} - \frac{w_{ij}^{(m-1)}}{x_i-x_j}\right),i,j=1,2,\dots,N;m=2,3,\dots,N-1
\tag{2.26}
\label{2.26}
$$
同样，$i=j$时的权值系数如下：
$$
\sum_{j=1}^N w_{ij}^{(m)}=0 或 w_{ii}^{(m)} = -\sum_{j=1,j\neq i}^N w_{ij}^{(m)}
\tag{2.27}
\label{2.27}
$$

## 3.广义积分求积法(GIQ)

广义积分求积(GIQ)是在与PDQ相同的概念下发展起来的。如果一个函数在整个区域内是光滑的，那么它可以用一个涉及整个区域内所有函数值的高次多项式来逼近。然后，通过对近似多项式进行积分，可以计算出函数在整个区域的一部分上的积分。因此，这种近似包含了整个区域中的所有函数值，而且精度很高，即使积分区域只包含两个点。显然，该方法的关键步骤是确定权重系数。下面将说明，GIQ的加权系数可以很容易地从PDQ中的一阶导数的加权系数中得到。

### 3.1  1维广义积分求积

假设$f(x)$在整个区域$[a,b]$上是连续的，且可以分解为$N-1$个区间，其网格点为$a=x_1,x_2,\dots,x_N=b$.由于$f(x)$在整个区域上是连续的，所以它可以用一个($N-1$)次多项式来逼近.特别地，当N个网格点的函数值已知时，$f(x)$可以用与所有网格点的函数值相关的拉格朗日插值多项式来逼近。因此，该逼近多项式在$[x_i,x_j]$上的积分可能涉及整环以外的函数值。作为一般情况，假设整个域的一部分上的$f(x)$的积分由整个域中的所有函数值与形式的线性组合来近似：
$$
\int_{x_i}^{x_j}f(x)\cdot dx = \sum_{k=1}^Nc_k^{ij} \cdot f(x_k)
\tag{3.1}
\label{3.1}
$$
与PDQ中的分析类似，作为$f(x)$的近似的($N-1$)次多项式构成了N维线性向量空间。因此，如果所有的基多项式都满足等式$(\ref{3.1})$，那么空间中的任何多项式都满足。考虑拉格朗日插值多项式$r_k(x),k=1,2,\dots,N$作为基函数，那么
$$
c_k^{ij} =\int_{x_i}^{x_j}r_k(x)\cdot dx
\tag{3.2}
\label{3.2}
$$
直接计算$c_k^{ij}$的表达是困难的，所以采用如下方法：
$$
f(x)=\frac{du(x)}{dx}
\tag{3.3}
\label{3.3}
$$
我们可以清楚地看到，如果$f(x)$是($N-1$)次多项式，则$u(x)$应该是$N$次多项式。假设$f(x)$由以下形式的($N-1$)次多项式逼近：
$$
f(x)=a_0+a_1x+a_2x^2+\dots+a_{N-1}x^{N-1}
\tag{3.4}\label{3.4}
$$
其中$a_0,a_1,a_{N-1}$为常数，对$(\ref{3.4})$从常数c到变量$x$作积分有：
$$
u(x)=\int_c^x f(t)dt+u(c)=F(x,c)+u(c)
\tag{3.5}\label{3.5}
$$
其中
$$
F(x,c)=x\left(a_0+\frac{a_1}{2}x+\dots+\frac{a_{N-1}}{N}x^{N-1}\right)-c\left(a_0+\frac{a_1}{2}c+\dots+\frac{a_{N-1}}{N}c^{N-1}\right)
\notag
$$
显然，$F(x,c)$取决于$N$个常数。证明了$F(x,c)$构成$N$维线性多项式向量空间。它的一组基多项式可选为:
$$
p_k(x)=(x-c)r_k(x),k=1,2,\dots,N
\tag{3.6}\label{3.6}
$$
相似于PDQ，可以得到：
$$
F_x(x_i,c)=\sum_{j=1}^N \underline{a}_{ij}F(x_j,c),i=1,2,\dots,N
\tag{3.7}\label{3.7}
$$
将$(\ref{3.6})$代入$(\ref{3.7})$中得到：
$$
\underline{a}_{ij}=\frac{x_i-c}{x_j-c}w_{ij}^{(1)},j\neq i
\tag{3.8a}\label{3.8a}
$$

$$
\underline{a}_{ii}=w_{ii}^{(1)}+\frac{1}{x_i-c},
\tag{3.8b}\label{3.8b}
$$

$$
\color{greenyellow}
p_k'(x)=r_k(x)+(x-c)r_k'(x)\\
\color{greenyellow}
\therefore~~(x_i-c)r_k'(x_i)=\underline{a}_{ij}(x_j-c) ~~~~~~~(i\neq  j)\\
\color{greenyellow}
\therefore~~\underline{a}_{ij}=\frac{x_i-c}{x_j-c}r_k'(x_i)=\frac{x_i-c}{x_j-c}w_{ij}^{(1)}\\
\color{greenyellow}
同样\\
\color{greenyellow}
p_k'(x_i)=1+(x_i-c)w_{ii}^{(1)} ~~~~~~~~~~~~(k=i)\\
\color{greenyellow}
\therefore \underline{a}_{ii}=w_{ii}^{(1)}+\frac{1}{x_i-c}
\notag
$$

从公式$(\ref{3.8a})-(\ref{3.8b})$中可以看出，不能选择c作为网格点的坐标，从公式10.6可以看出，当c=x时，积分区域为零。因此，不必近似$c=x$的积分。另一方面，公式$(\ref{3.8a})-(\ref{3.8b})$可以写成向量形式:
$$
\left\{F_x\right\}=\{\underline{A}\}\{F\}
\tag{3.9}
\label{3.9}
$$
其中
$$
\{F\} = \{F(x_1,c),F(x_2,c),\dots,F(x_N,c)\}^T\\
\{F_x\} = \{f\}=\{F_x(x_1,c),F_x(x_2,c),\dots,F_x(x_N,c)\}^T
\notag
$$
将$(\ref{3.4})$和$(\ref{3.6})$插入$(\ref{3.9})$中，令
$$
\{f^I\}=\{F\}=\{\int_c^{x_1}f(x)dx,\int_c^{x_2}f(x)dx,\dots,\int_c^{x_N}f(x)dx\}^T
\tag{3.10}
\label{3.10}
$$
可以得到：
$$
\{f\}=[\underline{A}]\{ f^I \}
\tag{3.11}
\label{3.11}
$$
令$[W^I]=[\underline{A}]^{-1}$，则有$\{f^I\}=W^I \{f\}$，进一步：
$$
\int_c^{x_i} f(x) dx =\sum_{k=1}^N w_{ik}^I f(x_k),~~~i=1,2,\dots,N
\tag{3.12}
\label{3.12}
$$

那么就得到
$$
c_k^{ij}=w_{jk}^I-w_{ik}^I
\tag{3.13}
\label{3.13}
$$
很明显GIQ中的权重系数是由PDQ中的一阶导数离散化得到的。
