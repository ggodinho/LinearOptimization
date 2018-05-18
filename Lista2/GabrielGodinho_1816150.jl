#Lista 2 - Implementando Simplex
#Gabriel Campos Godinho

maxIter = 1000 #Número máximo de iterações
tol = 1/10^(9)
problema = 0

function SimplexFase1(A,b,c,base,Nbase)
    println("Simplex Fase 1")
    m,n = size(A)
    A_w = -1*ones(m,n+1)
    A_w[:,1:n] = A
    b_w  = b
    c_w = zeros(n+1)
    c_w[n+1] = -1
    x = zeros(n+1)

    status = 0
    iter = 0

    while iter <= 1000
        iter = iter + 1
        println("Iter $(iter):")
        B = A_w[:, base]
        N = A_w[:, Nbase]
        cB = c_w[base]
        cN = c_w[Nbase]

        y = B'\cB
        cRed = (cN' - y'*N)'
        cRedMax,j = findmax(cRed)

        #Começando o processo iterativo
        xB = B\b
        x = zeros(Float64,n+1)
        x[base] = xB

        #Condição de convergência
        if cRedMax <= 0.0001
            status = 1
            println("x = $(x)")
            println("Base = $(base)")
            println("NBase = $(Nbase)\n")

            # reorder original matrix
            base_orig = [i for i in base if i!=(n+1)]
            Nbase_orig = [i for i in Nbase if i!=(n+1)]

            A1 = zeros(size(A))
            A1[:, 1:length(Nbase_orig)] = A_w[:,Nbase_orig]
            A1[:, (length(Nbase_orig)+1):n] = A_w[:, base_orig]
            c1 = zeros(n)
            c1[1:length(Nbase_orig)] = c_w[Nbase_orig]
            c1[(length(Nbase_orig)+1):end] = c_w[base_orig]
            return x,cB'*xB,status,base,Nbase
        end

        #direção extrema se ilimitado
        dB = -B\N[:,j]
        r = xB./(-dB)
        if all(r .< 0)
            d = zeros(n)
            d[Nbase] = 1
            d[base] = dB
            status = -1
            println("x = $(x)")
            println("Base = $(base)")
            println("NBase = $(Nbase)\n")
            return Inf, d, status
        end

        #Selecionando o menor valor não negativo para r
        r = xB./(-dB)
        k = length(r)
        for k=1:k
            if r[k] < 0
                r[k] = NaN
            end
        end
        i = indmin(r)

        println("x = $(x)")
        println("Base = $(base)")
        println("NBase = $(Nbase)\n")

        aux = base[i]
        base[i] = Nbase[j]
        Nbase[j] = aux
    end
end

function SimplexFase2(A,b,c,base,Nbase)
    m,n = size(A)
    status = 0
    println("Simplex Fase 2")
    iter = 0

    while iter <= maxIter
        iter = iter + 1
        println("Iter $(iter):")
        B = A[:,base]
        N = A[:,Nbase]
        cB = c[base]
        cN = c[Nbase]

        y = B'\cB
        cRed = (cN' - y'*N)'
        cRedMax,j = findmax(cRed)

        #Começando o processo iterativo
        xB = B\b
        x = zeros(Float64,n)
        x[base] = xB

        #Condição de convergência
        if cRedMax <= tol
            status = 1
            println("x = $(x)")
            println("Base = $(base)")
            println("NBase = $(Nbase)\n")
            return x,cB'*xB,status
        end

        #direção extrema se ilimitado
        dB = -B\N[:,j]
        r = xB./(-dB)

        if all(r .< 0)
            d = zeros(n)
            d[Nbase] = 1
            d[base] = dB
            status = -1
            println("x = $(x)")
            println("Base = $(base)")
            println("NBase = $(Nbase)\n")
            return Inf, d, status
        end

        #Selecionando o menor valor não negativo para r

        k = length(r)
        for k=1:k
            if r[k] < 0
                r[k] = NaN
            end
        end
        i = indmin(r)

        println("x = $(x)")
        println("Base = $(base)")
        println("NBase = $(Nbase)\n")

        aux = base[i]
        base[i] = Nbase[j]
        Nbase[j] = aux
    end
end

function printSolution(x, z, status)
    println("\nSolução obtida:")
    if status == -1
        println("Problema ilimitado")
        println(("d = $(z)"))
    else
        println("x = $(x)")
        println("z = $(z)")
    end
    println("status = $(status)")
end

for problema = 1:3

    if  problema == 1 #Problema 1
        A = [2 1 1 0; 1 2 0 1]
        b = [4;4]
        c = [4; 3; 0; 0]
    elseif problema == 2 #Problema 2
        A = [0.5 -1 1 0; -4 1 0 1]
        b = [0.5;1]
        c = [1; 1; 0; 0]
    elseif problema == 3 #Problema 3
        A = [2 1 1 0 0; 1 2 0 1 0; -1 -1 0 0 1]
        b = [4;4;-1]
        c = [4; 3; 0; 0; 0]
    end

    println("\n========================== Problema $(problema) ==============================")
    m,n = size(A)
    base = [i for i = (n-m+1):(n)]
    Nbase = [i for i = 1:(n-m) if !(i in base)]

    if minimum(b) < 0
        #Preparando nova base
        p = indmin(b)
        base[p] = n+1
        Nbase = vcat(Nbase, n-m + p)

        #Resolvendo Simplex fase 1
        x1,z1,status1,base1,Nbase1 = SimplexFase1(A,b,c,base,Nbase)

        i = indmax(Nbase1)
        Nbase1 = deleteat!(Nbase1,i)

        #Resolvendo o Simplex fase 2
        x,z,status = SimplexFase2(A,b,c,base1,Nbase1)
        printSolution(x,z,status)
    else
        x,z,status = SimplexFase2(A,b,c,base,Nbase)
        printSolution(x,z,status)
        println("")
    end
end
