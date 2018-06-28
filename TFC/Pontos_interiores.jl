#Algoritmo de pontos interiores - TFC Programação Linear

function PInteriores(A,b,c,x,y,s)
    #Inicializando
    solution = 0
    m,n = size(A)
    k = 0
    erro = 1/10^5
    iter_max = 100
    x_1 = zeros(iter_max+1)
    x_2 = zeros(iter_max+1)

    for k = 0:iter_max
        #Teste de otimalidade
        viab_primal = maximum(A*x - b) < erro
        viab_dual = maximum(A'y - s - c) < erro

        if abs(s'*x) < erro
            solution = 1
            x_1 = x_1[1:k]
            x_2 = x_2[1:k]
            return solution,x,y,s,k,viab_primal,viab_dual,x_1,x_2
            break
        end

        #Metodo de Newton
        p = 0.9
        u = p*(s'x)/n #Cálculo de u (mi)
        X = diagm(x,0)
        S = diagm(s,0)

        d = [A zeros(m,m) zeros(m,n); zeros(n,n) A' diagm(ones(n),0); S zeros(n,m) X]\
            -[A*x - b; A'y - s - c; X*S*ones(n) - u*ones(n)]

        d_x = d[1:n]
        d_y = d[n+1:n+m]
        d_s = d[n+m+1:length(d)]

        #Encontrando o menor d_x e d_s para cálculo de beta
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

        #Cálculo de beta
        beta_p = min(1,(0.9*(-x[ind]/-aux_dx)))
        beta_d = min(1,(0.9*(-s[ind_s]/-aux_ds)))

        x_1[k+1] = x[1]
        x_2[k+1] = x[2]

        #Atualização de x, y e s com os passos calculados
        x = x + beta_p*d_x
        y = y + beta_d*d_y
        s = s + beta_d*d_s
    end
end

#Problema da produção
A = [2 1 1 0; 1 2 0 1]
b = [4;4]
c = [4; 3; 0; 0]

#Escolha de variáveis iniciais
x = [1; 1; 5; 5]
y = [5; 2]
s = [5;2;5;1]

solution,x,y,s,k,v_primal,v_dual,x1_path,x2_path = @time PInteriores(A,b,c,x,y,s)
z = c'x

#Escrevendo no arquivo texto o caminho percorrido por x
open("C:/Users/Gabriel_2/Dropbox/Mestrado/2018-01-Programação Linear/GithubRep/TFC/path.txt", "w") do f
    write(f, "$x1_path, \n", "$x2_path")
end
