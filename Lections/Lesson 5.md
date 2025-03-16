Храним только ненулевые элементы матриц, можем работать с очень большими размерностями.
$\tau = 1/(\lambda_{min} + \lambda_{max})$
$q = (\varkappa(A)- 1)/(\varkappa(A) + 1)$
$\varkappa(n) \sim c*N^2$

1. Накопление ошибки
	$x_k -> x_{k+1}$
2. $Г_k \neq$
![[IMG_9655.jpg]]
![[IMG_9656.jpg]]
![[IMG_9657.jpg]] 
![[IMG_9658.jpg]]
![[IMG_9659.jpg]]
![[IMG_9662.jpg]]
![[IMG_9663.jpg]]
![[IMG_9665.jpg]]
![[IMG_9666.jpg]]
![[IMG_9668.jpg]]
![[IMG_9671.jpg]]
![[IMG_9673.jpg]]
![[IMG_9674.jpg]]
![[IMG_9676.jpg]]
![[IMG_9676 (2).jpg]]
![[IMG_9681.jpg]]
![[IMG_9683.jpg]]
## Выбор предобуславливателя
1. Якоби
	M = $diag(A)$
2. Гаусса-Зейделя
	Нижний или верхний треугольник, включая диагональ
	$Ax=b$
	$Lx + Dx + Rx = b$
	$Lx_{k+1} + Dx_{k+1} + Rx_k = b$
	$(L+D)(x_{k+1}-x_k) + \underbrace{Lx_k + Dx_k + Rx_k}_A = b$
	$(L+D)\frac{x_{k+1}-x_k)}{1} + Ax_k = b$
	$Lx_{k} + Dx_{k+1} + Rx_{k+1} = b$
	$\underbrace{(D+R)}_M\frac{x_{k+1}-x_k)}{1} + Ax_k = b$
	$\tau = 1$
	Но эта матрица несимметрична. Поэтому сделаем следующим образом:
	$x_k \rightarrow x_{k+1}$
	$Lx_{k+1/2} + Dx_{k+1/2} + Rx_k = b$
	$Lx_{k+1/2} + Dx_{k+1} + Rx_{k+1} = b$
	-- Симметричный метод Зейделя
	$x_{k + 1/2} = (L + D)^{-1}(b - Rx_k)$
	$x_{k + 1} = (D+R)^{-1}(b - L(L + D)^{-1})(b - Rx_k)=$
	$=(D+R)^{-1}(b - L(L + D)^{-1})(b - Rx_k)$
	$(D+R)(x_{k+1} - x_{k + 1/2}) + Ax_{k+1/2} = b$
	$\Rightarrow x_{k+1} = x_{k+1/2}(D+R)^{-1}(Ax_{k+1/2}-b)$
	$(L+D)(x_{k+1} - x_{k + 1/2}) + Ax_k=b$
	$x_{k+1/2} = x_k - (L+D)^{-1} \frac{r_k}{Ax-b}$
	$x_{k+1}=x_k - (L+D)^{-1}r_k - (D+R)^{-1}(A(x_k-(L+D)^{-1}r_k) - b)$
	$=x_{k}-((D+L)^{-1}+(D+R)^{-1} - (D+R)^{-1}A(L+D))r_k$
	$=x_{k}-(D+R)^{-1}((D+R)(D+L)^{-1}+E-A(L+D))r_k$
	$=x_{k}-(D+R)^{-1}(D+R+L+D-\underbrace{A}_{D+R+L})(L+D)^{-1}r_k$
	$=x_{k}-(D+R)^{-1}D(L+D)^{-1}r_k$
	$x_{k+1}-x_{k}=(D+R)^{-1}D(L+D)^{-1}(Ax_k-b)$
	$\underbrace{(D+L)D^{-1}(D+R)}_{M}(x_{k+1}-x_k)+Ax_k=b$
	$M^T=M$
	$M$ - предобуславливатель Зейделя
	$Mx=y$
	a) $(D+L)x_1=y$
	b)$x_2 = Dx_1$
	c) $(D+R)x_3=x_2$
3. Метод верхней релаксации
	$A = L + D + R$ - не единственный способ представить матрицу в виде суммы обратимых
	$\omega A = \omega(L+D+R)=(D+\omega L) + (\omega R - (1 - \omega)D)$
	$\omega \in (0,2)$
	$\omega = 1$ - Зейдель
![[IMG_9684.jpg]]
![[IMG_9685.jpg]]
![[IMG_9687.jpg]]
![[IMG_9691.jpg]]
![[IMG_9697.jpg]]
![[IMG_9698.jpg]]
![[IMG_9700.jpg]]

## Неполное разложение Холецкого
$A = R^tDR$, D - diaf(+-1), R - верхнетреуг
$a_{ij} = \sum_{k=1}^{n}(R^t)_{ik}d_{kk}r_{kj} = \sum_{k=1}^{n}r_{kj}d_{kk}r_{kj}$
$a_{ij} = \sum_{k=1}^{min(i,j)}r_{ki}d_{kk}r_{kj}$
$a_{ii} = \sum_{k=1}^{i}r_{ki}d_{kk}r_{ki} = \sum_{k=1}^{i-1}r_{ki}^2d_{kk} + d_{kk}r_{ii}^2$
$r_{ii} = \sqrt{|a_{ii} - \sum_{k=1}^{i-1}r_{ki}^2d_{kk}|}$
$d_{ii} = sign(|a_{ii} - \sum_{k=1}^{i-1}r_{ki}^2d_{kk}|)$
$a_{ij} = \sum_{k=1}^{i}r_{ki}d_{kk}r_{kj} = \sum_{k=1}^{i-1}r_{ki}d_{kk}r_{kj} + r_{ii}d_{ii}r_{ij}$
$\Rightarrow r_{ij} = (a_{ij} - \sum_{k=1}^{i-1}r_{ki}d_{kk}r_{kj}) / (r_{ii}d_{ii})$
