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
	$n = 2:W_2^1(\Omega) \subset L_{\infty}(\Omega)$
	$n = 3:W_2^1(\Omega) \subset L_{n}(\Omega)$
8. $W_p^l (\Omega) = \{u \in L_p(\Omega): \exists \frac{\partial u}{\partial x_i}\}$
9. Неравенства Фридрихса
	$\exists c_1 = c_1(\Omega): \forall u \in W_2^1(\Omega), u|_{\partial \Omega} = 0$
	$||u||_{L_2(\Omega)} \leq c_1 ||\nabla u||_{L_2(\Omega)}$
	$||u||^2_{W_2^1(\Omega)}$
	![[IMG_9903.jpg]]
10. Неравенства Пуанкаре
	$\exists c_2 = c_2(\Omega): \forall u \in W_2^1(\Omega) \quad (\nabla u, n)|_{\partial \Omega} = 0$
	$\int_{\Omega}udx = 1 \quad (u, n) = 1$
	$||u||_{L_2(\Omega)} \leq C_2 ||\nabla u||_{L_2(\Omega)}$
	$W_2^1(\Omega) / \Real = \{u \in W_0^1(\Omega):(\)$
	![[IMG_9904.jpg]]

$$
\begin {cases}
-div(A \nabla u) = f \\
u|_{\partial \Omega} = 0
\end {cases} \\
\begin {cases}
-div(A \nabla u) = f \\
(A \nabla u, n)|_{\partial \Omega} = 0
\end {cases} \\
$$

![[IMG_9905.jpg]]

![[IMG_9907.jpg]]

Получили обобщённое решение уравнения аля Лапласса

![[IMG_9908.jpg]]

Определили решение, когда $f \in L_2(\Omega), A = (a_ij) \in L_{\infty}, \Omega$ - измерима
1. Корректность
	1. Совпадение настоящ. и обобщ. п. в., если настоящ. уществует
	2. $\f \in L_2, A \in L^{\infty}, \partial\Omega - C^1$ - тогда совп почти всюду
$(A \nabla u, \nabla v)_{L_2(\Omega)} = (f, v)_{L_2(\Omega)}$
$\lambda_{min}(x,x) \leq (Ax,x) \leq \lambda_{max}(x,x)$
$\lambda_{min}(\nabla x,\nabla x) \leq (A\nabla x,\nabla x) \leq \lambda_{max}(\nabla x,\nabla x)$
$(\nabla u,\nabla v) = ||u||_{W_2^1(\Omega)}$