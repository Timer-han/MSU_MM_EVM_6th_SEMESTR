Рассмотрим решение диффура
$$
\Delta u = u_{x_1x_1} + u_{x_2x_2}
$$
$$
-\Delta u = f
$$
$$
u|_{\partial \Omega} = f
$$
$$
div(A, \nabla u)
$$
$$
(a_{ij} \in L_{\infty(\Omega)})
$$
Дирихле 
$$
\begin {cases}
-div(A, \nabla u) = f \\
u|_{\partial \Omega} = g
\end {cases}
$$

![[IMG_9895.jpg]]

![[IMG_9896.jpg]]

![[IMG_9897.jpg]]



Гаус-Остроградский, Грин

![[IMG_9899.jpg]]

![[IMG_9900.jpg]]

![[IMG_9901.jpg]]

![[IMG_9902.jpg]]

1. Корректно определена
	1. Единственность
	2. Если существует простая, то обобщённая производная совпадает с настоящей
2. Свойства
	$W_2^1(\Omega) = \{u \in L_2(\omega): \exists \frac{\partial u}{\partial x_i} \in L_2(\Omega)\}$
	$(u, v)_{W_2^1(\Omega)} = (u, v)_{L_2(\Omega)} + (\nabla u, \nabla v)_{L_2(\Omega)}$
3. $||u||_{W_2^1(\Omega)} = (u, u)_{W_2^1(\Omega)}^{\frac{1}{2}}$
4. Полное пространство
5. $\int_{\partial \Omega}f(a, n)dS = \int_{\Omega}fdiv(a)dx + \int_{\Omega}(a, \nabla f)dx$
	$\int_{\partial \Omega}f\cdot u \cdot dS = \int_{\Omega}f \cdot \frac{\partial u}{\partial x_i}dx + \int_{\Omega}u \frac{\partial f}{\partial x_i}dx$
6. Теорема о сладах
	$\frac{\partial u}{\partial n} = (\nabla u, \overrightarrow n)$
7. Соболев. Теорема о вложении
	$W_2^1(\Omega) \subset L_2(\Omega)$
	$W_2^1(\Omega) \subset L_p(\Omega) \forall p \geq 2 ?$
	$W_2^1(\Omega) \subset L_p(\Omega)$
