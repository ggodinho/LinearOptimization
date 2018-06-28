#Problema da produção
A = [2 1 1 0; 1 2 0 1]
b = [4;4]
c = [4; 3; 0; 0]

#Inicializando
solution = 0
m,n = size(A)
k = 0
x = [10; 20; 5; 5]
y = [5; 5]
s = [5;5;5;5]
erro = 1/10^5
iter_max = 100
x_1 = zeros(iter_max+1)
x_2 = zeros(iter_max+1)
beta_p = 0.0
beta_d = 0.0
d_x = 0.0
viab_primal = A*x - b == 0
viab_dual = A'y - s - c

for k = 0:iter_max
        #Teste de otimalidade
    viab_primal = maximum(A*x - b) < erro
    viab_dual = maximum(A'y - s - c) < erro

    if abs(s'*x) < erro
        solution = 1
        break
    end

    #Metodo de Newton
    rho = 0.9
    mi = rho*(s'x)/n
    X = diagm(x,0)
    S = diagm(s,0)
    id = [1]

    d = [A zeros(m,m) zeros(m,n); zeros(n,n) A' diagm(ones(n),0); S zeros(n,m) X]\
        -[A*x - b; A'y - s - c; X*S*ones(n) - mi*ones(n)]

    d_x = d[1:n]
    d_y = d[n+1:n+m]
    d_s = d[n+m+1:length(d)]

    aux_dx = zeros(n)
    for i = 1:n
        if d_x[i]>=0
            aux_dx[i]= 0
        else
            aux_dx[i] = -d_x[i]
        end
    end
    ind = indmax(aux_dx)
    aux_dx = maximum(aux_dx)

    aux_ds = zeros(n)
    for i = 1:n
        if d_s[i]>=0
            aux_ds[i]= 0
        else
            aux_ds[i] = -d_s[i]
        end
    end
    ind_s = indmax(aux_ds)
    aux_ds = maximum(aux_ds)

    #beta_p = 0.9*min(1.0,maximum(-x/minimum(min(0.0,d_x))))
    beta_p = min(1,(0.9*(-x[ind]/-aux_dx)))
    beta_d = min(1,(0.9*(-s[ind_s]/-aux_ds)))

    x_1[k+1] = x[1]
    x_2[k+1] = x[2]

    x = x + beta_p*d_x
    y = y + beta_d*d_y
    s = s + beta_d*d_s
end
