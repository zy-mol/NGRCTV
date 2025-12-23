# NGRCTV
## Optimization Model

$$\begin{align*}
\min_{\mathcal{U},\mathbf{V},\mathcal{E},\mathcal{S}} & \tau_1\sum_{k=1}^R\| \mathcal{G}_1(\cdot,\cdot,k)^\top\|_{2,1}+\tau_2\sum_{k=1}^R\| \mathcal{G}_2(\cdot,\cdot,k)\|_{2,1}+\beta\|\mathcal{E}\|_F^2+\lambda\|\mathcal{S}\|_1,\\
& \textrm{s.t. } \nabla_n(\mathcal{U})=\mathcal{G}_n,\, n=1,2,\\
& \mathcal{Y}=\mathcal{U}\times_3\mathbf{V}^\top+\mathcal{E}+\mathcal{S}, \ \mathbf{V}^\top\mathbf{V}=\mathbf{I}.
\end{align*}$$
