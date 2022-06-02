using Revise

envnames = ["period", "muperiod", "eigen"]

# make directories to store log files
_ckmkdir(f) = isdir(f) == false && mkdir(f)
_ckmkdir("log")
for envname in envnames
    _ckmkdir("log/"*envname)
end

#--------------- for 4 seed numbers ---------------------# 
seednum = 1234
include("scripts/period.jl")
seednum = 2345
include("scripts/period.jl")
seednum = 3456
include("scripts/period.jl")
seednum = 4567
include("scripts/period.jl")

seednum = 1234
include("scripts/period_2.jl")
seednum = 2345
include("scripts/period_2.jl")
seednum = 3456
include("scripts/period_2.jl")
seednum = 4567
include("scripts/period_2.jl")

seednum = 1234
ConstN = 1.
include("scripts/eigen.jl")
seednum = 1234
ConstN = 5.
include("scripts/eigen.jl")
seednum = 1234
ConstN = 10.
include("scripts/eigen.jl")
seednum = 1234
ConstN = 30.
include("scripts/eigen.jl")

#--------------- compute the last results ---------------# 

for envname in envnames
    # reading data
    d1 = read("log/"*envname*"/"*envname*".jlso", JLSOFile)
    d2 = read("log/"*envname*"/"*envname*"_1.jlso", JLSOFile)
    d3 = read("log/"*envname*"/"*envname*"_2.jlso", JLSOFile)
    d4 = read("log/"*envname*"/"*envname*"_3.jlso", JLSOFile)
    if envname == "period"
        estimation1 = sum(d1[:logs][:algstate][:est])
        estimation2 = sum(d2[:logs][:algstate][:est])
        estimation3 = sum(d3[:logs][:algstate][:est])
        estimation4 = sum(d4[:logs][:algstate][:est])
        println("Estimation for $envname -> 1: $estimation1 ,  
        2: $estimation1 ,  3: $estimation3 ,  4: $estimation4")
    elseif envname == "muperiod"
        estimation1 = sum(d1[:logs][:algstate][:est])
        estimation2 = sum(d2[:logs][:algstate][:est])
        estimation3 = sum(d3[:logs][:algstate][:est])
        estimation4 = sum(d4[:logs][:algstate][:est])
        println("Estimation for $envname -> 1: $estimation1 ,  
        2: $estimation1 ,  3: $estimation3 ,  4: $estimation4")
    elseif envname == "eigen"
        matrix1 = Matrix(d1[:buffers].lastoutput)
        matrix2 = Matrix(d2[:buffers].lastoutput)
        matrix3 = Matrix(d3[:buffers].lastoutput)
        matrix4 = Matrix(d4[:buffers].lastoutput)
        eig1 = eigen(matrix1).values
        eig2 = eigen(matrix2).values
        eig3 = eigen(matrix3).values
        eig4 = eigen(matrix4).values
        println("Estimation of eigenvalues for $envname -> \n\n1: $eig1 \n\n
        2: $eig2 \n\n  3: $eig3 \n\n 4: $eig4")
    end
end
