
# first run (strg+enter) the cell at the very bottom of this notebook to define all functions

# the output files are
# a) the coordinates of the final configuration (x1,y1,z1,x2,y2,z2,x3,...)
# b) the final configuration's radius of gyration
# c) the final adjacency matrix. '1' means linked, '2' means not linked

#here we define all parameters
#default values correspond to the chain of equally-sized spheres described in the main text

N=10              #number of spheres
σ=0.0             #variation in sphere radii. radius = 0.5* (1+σ*RandomNumber∈[0,1]) for each sphere independently
K=50              #spring constant for harmonic potential. The higher, the more accurately are the constraints fulfilled, but the smaller the numerical time step must be
dt=0.8*1e-2       #numerical time step 
ΔT=0.15           #time interval after which shrinking the distance between chosen spheres is checked
A=5               #"bending stiffness" of initial chain. The larger, the straighter.
NumberOfSims=30; #how many independent simulations are run for this parameter set

#start simulations. output files are written at the bottom of this cell.
for sim in 1:NumberOfSims

tic()

r=0.5*(ones(N)+σ*(rand(N)-0.5)*2)  #sphere radii for THIS simulation
d=r[1:N-1]+r[2:N]                  #corresponding bond lengths
R0=SetIC(N,A,d);                   #draw initial condition from Boltzmann distribution
    
M=zeros(N,N)                       #generate adjacency matrix
M[1,1]=M[1,2]=M[N,N]=M[N,N-1]=2    #set (next) diagonal to 1
for i=2:N-1
  M[i,i-1]=M[i,i]=M[i,i+1]=2
end


#First Bond
pair=FindPair(M)                   #Choose untried entry (=0) in M
p=pair[1];q=pair[2]

                                   #first one manually for programming syntax reasons ('bonds')
R0=Simulate(R0,dt,ΔT,d,p,q,K,[],r) #run simulations in which p and q are pulled together
M[p,q]=M[q,p]=1                    #assume first one worked.
bonds=[p q]                        #keeps track of the successfully formed bonds


for i in 2:(N^2-N-2*(N-1))/2       #loop over number of remaining zeros in M
    pair=FindPair(M)               #Choose untried entry (=0) in M
    p=pair[1];q=pair[2]
    Rstart=R0;                     #remember state before this cycle
    R0=Simulate(R0,dt,ΔT,d,p,q,K,bonds,r)                      #run simulations in which p and q are pulled together
    if norm(R0[(3p-2):(3p)]-R0[(3q-2):(3q)]) < r[p]+r[q]       #check distance between p and q after simulation.
        M[p,q]=M[q,p]=1            #in case of success matrix entry to '1' and add new bond
        bonds=[bonds; [p q]]
    else
        M[p,q]=M[q,p]=2            #in case of no success set matrix entry to '2' and restore configuration before this cycle
        R0=Rstart;
    end

        

    if mod(i,100)==0               #every 100th cycle, we check if we can already exclude some untried bonds
            
        # check whether there are more than 8N ones in M, as this would mean we've achieved the 'closest ball packing' possible in 3D
        MatCount=0
        for element in M
            if element==1
               MatCount+=1
            end
        end
        if MatCount>(8*N)-1
            for element in M
                if element==0
                   element=2
                end
            end
        end

        # check whether any bead is bonded to 12 other beads, as this is the maximum number of same-sized spheres a sphere can touch.
        for row in 1:N
            LineCount=0
            for element in M[row,:]
                if element==1
                    LineCount+=1
                end
            end
            if LineCount>11
                for element in M[row,:]
                    if element==0
                        element=2
                    end
                end
            end
        end
    end
end

#write results to files
writedlm("results/AdjacencyMatrix_N$(N)_$(sim).txt",M)
writedlm("results/FinalCoordinates_N$(N)_$(sim).txt",R0)
writedlm("results/RadiusOfGyration_N$(N)_$(sim).txt",Rgyration(R0,1,N))

toc()
    
end

# total force vector
function f(R,d,p,q,K,bonds,NB,r)
    
    return ForceConnect(R,p,q)+K*(ForceExvol(R,r,NB)+ForceBonds(R,r,bonds)+ForceStretch(R,d))
    
end



#find an untried pair of spheres
function FindPair(M)

N=size(M,1)
ValidIndex=[]
for i=1:N^2
  if M[i]==0
     append!(ValidIndex,i)
  end
end

x=rand(ValidIndex)
i=Integer(floor((x-1)/N))+1
j=x-N*(i-1)

return[i,j]
    
end



#forces that hold spheres together which were successfully connected earlier in the simulation
function ForceBonds(R,r,bonds)

BondNumber=Integer(size(bonds,1))
N=Integer(length(R)/3)
Force=zeros(3N)
 
for i=1:BondNumber
   p=bonds[i,1]
   q=bonds[i,2]
   vec=R[(3p-2):(3p)]-R[(3q-2):(3q)]
   dist=norm(vec)
   WW=(dist-(r[p]+r[q]))*vec/dist
   Force[(3p-2):(3p)]-=WW
   Force[(3q-2):(3q)]+=WW
end

return Force

end



# force that pulls the two chosen spheres together
function ForceConnect(R,p,q)
    
Force=zeros(3N)
univec=R[(3p-2):(3p)]-R[(3q-2):(3q)]
r=norm(univec)
univec/=norm(univec)

Force[(3p-2):(3p)]-=univec
Force[(3q-2):(3q)]+=univec
    
return Force

end



#excluded volume forces
function ForceExvol(R,r,NB)

N=Integer(length(R)/3)
    
Force=zeros(3N)

for i in 1:N
    for j in NB[i]
      if abs(i-j)>1
        rij=R[(3i-2):(3i)]-R[(3j-2):(3j)]
        dij=norm(rij)
        if dij<r[i]+r[j]
          WW=(r[i]+r[j]-dij)*rij./dij
          Force[(3i-2):(3i)]+=WW
          Force[(3j-2):(3j)]-=WW
        end
        deleteat!(NB[j],findin(NB[j],i))
      end
    end 
end

return Force;
    
end



#stretching forces, holding together spheres that are connected along the backbone
function ForceStretch(R,d)

N=Integer(length(R)/3)
Force=zeros(3N)
 
for i=1:N-1
   vec=R[(3i+1):(3i+3)]-R[(3i-2):(3i)]
   dist=norm(vec)
   WW=(dist-d[i])*vec/dist
   Force[(3i-2):(3i)]+=WW
   Force[(3*i+1):(3*i+3)]-=WW
end

return Force

end



#updating the neighbour lists which, for each sphere, keep track of spheres spatially close to it
function NBList(R,threshold)

N=Integer(length(R)/3)

NB = Array{Any}(N)
for index in eachindex(NB)
  NB[index]=Int64[]
end

for i in 1:N
  for j in 1:N
    if i != j
      if norm(R[(3i-2):(3i)]-R[(3j-2):(3j)])<threshold
        append!(NB[i],j)
      end
    end
  end
end   

return NB
    
end



#check for overlaps during generating an initial chain
function Overlap(R,x)

N=Integer(length(R)/3)
    
if x == 0
   return 0
end

for i in 3:N
    for j in 1:i-2
        if norm(R[(3i-2):(3i)]-R[(3j-2):(3j)])<2x
            return 1
        end
    end 
end

return 0
    
end



#radius of gyration for subchain starting at Bead1, ending at Bead2
function Rgyration(Positions,Bead1,Bead2)

NumberOfTimes=size(Positions,2)
NumberOfBeads=Bead2-Bead1+1
Positions=Positions[1:3*NumberOfBeads,:]
TotalSum=0;

for iTime in 1:NumberOfTimes     #Summing up 1/(N-1)*sum_Beads (Vec_i-COM)^2, Vec_i=position vector of bead i
   COM=[mean(Positions[1:3:end-2,iTime]),mean(Positions[2:3:end-1,iTime]),mean(Positions[3:3:end,iTime])]
   Sum=0
   for iBead in 1:NumberOfBeads
      Vec=Positions[(3*iBead-2):(3*iBead),iTime]-COM
      Sum+=Vec'*Vec
   end
   TotalSum+=Sum/NumberOfBeads
end
    
return sqrt(TotalSum/NumberOfTimes)         #Return average
    
end



#standard explicit Runge-Kutta method of 4th order
function RK4(R,dt,d,p,q,K,bonds,NB,r)
    
    k1=f(R,d,p,q,K,bonds,NB,r)
    k2=f(R+0.5*k1*dt,d,p,q,K,bonds,NB,r)
    k3=f(R+0.5*k2*dt,d,p,q,K,bonds,NB,r)
    k4=f(R+k3*dt,d,p,q,K,bonds,NB,r)
    
    return dt*(k1/6+k2/3+k3/3+k4/6)
    
end



#rotate a 3D vector around an axis by an angle α
function rotate(Vector,RotationAxis,α)

RotationMatrix=zeros(3,3);
n1=RotationAxis[1];n2=RotationAxis[2];n3=RotationAxis[3];

a=sin(α);b=cos(α);

RotationMatrix[1,1]=n1^2*(1-b)+b;
RotationMatrix[1,2]=n1*n2*(1-b)-n3*a;
RotationMatrix[1,3]=n1*n3*(1-b)+n2*a;
RotationMatrix[2,1]=n2*n1*(1-b)+n3*a;
RotationMatrix[2,2]=n2^2*(1-b)+b;
RotationMatrix[2,3]=n2*n3*(1-b)-n1*a;
RotationMatrix[3,1]=n3*n1*(1-b)-n2*a;
RotationMatrix[3,2]=n3*n2*(1-b)+n1*a;
RotationMatrix[3,3]=n3^2*(1-b)+b;

return RotationMatrix*Vector
    
end



#generate initial chain
function SetIC(N,A,d)
    
R=[0.;0.;0.]
limit=9*N
    
while length(R)<3N
        
R=[0.;0.;0.]
T=[0.;0.;d[1]]
fac=exp(A)-exp(-A)
counter=1
Iteration=0
limit+=N

while counter<N && Iteration<limit
    
  num=rand(1)[1]
  cosine=1/A*log(num*fac+exp(-A))
  Vec1=cross(T,T+rand(3));Vec1/=norm(Vec1)
  Vec2=rotate(Vec1,T,2*pi*rand(1)[1])
  T=rotate(T,Vec2,acos(cosine));T=T/norm(T)
  NewBead=R[end-2:end]+T*d[counter]

  if Overlap(R,0.5) == 0
    append!(R,NewBead)
    counter+=1
  end

  Iteration+=1

end
    
end

R[1:3:end-2]-=mean(R[1:3:end-2])
R[2:3:end-1]-=mean(R[2:3:end-1])
R[3:3:end]-=mean(R[3:3:end]);
    
return R

end



#simulating a cycle until the chosen spheres meet, or not move any more (see main text)
function Simulate(R,dt,ΔT,d,p,q,K,bonds,r)

N=Integer(length(R)/3)
NextCheck=ΔT
Time=0
vec=R[(3p-2):(3p)]-R[(3q-2):(3q)]
dpq=norm(vec)
LastR=dpq
NB=NBList(R,1.2)
maxdistance=0
    
while dpq>r[p]+r[q]
        
    #integrate
    inc=RK4(R,dt,d,p,q,K,bonds,1*NB,r)
    R+=inc
    maxdistance+=abs(maximum(inc))*sqrt(3)
    if maxdistance>0.2
      maxdistance=0
      NB=NBList(R,2*maximum(r)+0.2)
    end

    #check whether distance has decreased sufficiently
    if Time>NextCheck
        NextCheck+=ΔT
        NewR=norm(R[(3p-2):(3p)]-R[(3q-2):(3q)])
        if NewR>LastR-0.3*2*ΔT/(N/2)
           break
        end
        LastR=NewR
    end         
        
    #update time
    Time+=dt

    vec=R[(3p-2):(3p)]-R[(3q-2):(3q)]
    dpq=norm(vec)
end
    
    
return R
    
end
