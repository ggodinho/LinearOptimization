# Simplex Versão 2
println("")
println("")

println("===== Modelo de Otimização =====")

#Problema que inventei
A = [1 -1 1 0 0; 1.5 1 0 1 0; 0 1 0 0 1]
b = [-1, 3, 2]
c = [3 0 0 0 0]

# # Problema 1
# A = [2 1 1 0; 1 2 0 1]
# b = [4, 4]
# c = [4 3 0 0]

#Problema 2
# A = [0.5 -1 1 0; -4 1 0 1]
# b = [0.5; 1.0]
# c = [1 1 0 0]

#Problema caderno
# A = [3 3 1 0 0; 1 2 0 1 0; -1 -1 0 0 1]
# b = [4, 4, -1]
# c = [4 3 0 0 0]

SimplexCompleto(A,b,c)

function SimplexCompleto(A,b,c)
    m,n2 = size(A)
    n = n2 - m
    base = Array((n+1):(n + m))
    nbase = Array((1):(n))

    if all(b .> 0)
        println("")
        println("===== SIMPLEX - Fase 2 =====")
        println("")
        x,z,status = Simplex(A,b,c,base,nbase)
    else
        println("")
        println("===== SIMPLEX - Fase 1 =====")
        println("")

        i = indmin(b)
        base[i] = n + m + 1
        nbase = vcat(nbase, n + i)

        c_1 = zeros(n + m + 1)
        c_1[n + m + 1] = -1
        A_1 = hcat(A, -ones(m))

        aux1,aux2,aux3,base_2,nbase_2 = Simplex(A_1,b,c_1,base,nbase)


        println("")
        println("===== SIMPLEX - Fase 2 =====")
        println("")

        i = indmax(nbase_2)
        nbase_2 = deleteat!(nbase_2,i)

        x,z,status = Simplex(A,b,c,base_2,nbase_2)

    end
end


function Simplex(A,b,c,base,nbase)
    m = size(A)[1]
    n = size(A)[2]

    iter = 1
    while iter <= 20
        println("Iteração: ", iter)
        B, N = A[:,base], A[:,nbase]
        cb, cn = c[base], c[nbase]

        println("Decisão Atual: ", base)


        c_red = cn' - cb' * (B \ N)
        j = indmax(c_red)

        xb = B \ b
        x = zeros(n)
        x[base] = xb

        if all(c_red .<= 0)
            z = c[base]' * x[base]
            println("Custo Total = ", z)
            println("Variáveis Básicas / Valores = ", base, "/", x[base])
            status=1
            return x,z,status,base,nbase
        end

        y = -B \ N[:,j]
        r = x[base] ./ y

        d = zeros(n)
        d[nbase[j]] = 1
        d[base] = y


        if all(r .> 0)
            println("Problema Irrestrito.")
            println("Custo Total = inf")
            println("Direção Extrema = ", d)
            z = "inf"
            status = -1
            return x,z,status
        end

        aux2 = 0
        valor,i = findmax(r)
        while valor > 0 & aux2 <= m
            aux2 += 1
            # println("Valor de r = ", r[i])
            if valor >= 0.
                r[i] = NaN
            end
            valor,i = findmax(r)
        end

        aux = base[i]
        base[i] = nbase[j]
        nbase[j] = aux
        iter += 1
        println("")
    end

end
