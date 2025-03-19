#  Function SIS let us obtain  values for:
#  S_n,I_n
#  T is the steps number.

using Random, Distributions

function SDE_SIS(T::Int64, Time::Int64, h::Float64, mean1::Float64, desvest::Any, gamma::Float64, beta::Any, I0::Float64)
    #SIS(T,Time,h,mean,desvest,gamma,beta,I0)

    #initial disease rate
    bt = beta[1]
    desv = desvest[1]
    S1 = zeros(T + 1, 2)#auxiliar intermediate steps
    #create_vectors of returns
    S = zeros(Time)
    I = zeros(Time)

    #Initialization of I_0 S_0 
    S1[1, 2] = I0[1] #I_0
    S1[1, 1] = 1.0 - I0[1] #S_0
    # SDE SIS FIRST EVOLUTION
    for j in 2:(T+1)
        z = rand(Normal(mean1, 1))
        #evaluate Infected next  step j - rate infection with additive noise evaluating with step j-1 data
        S1[j, 2] = S1[j-1, 2] + (bt * h + desv * z * sqrt(h)) * (1.0 - S1[j-1, 2]) * S1[j-1, 2]  - gamma * S1[j-1, 2] * h#Infected
        S1[j, 1] = 1.0 - S1[j, 2]
    end
    #saving the first values -  final subdivision
    S[1] = S1[1, 1]
    I[1] = S1[1, 2]

        for j in 2:(T+1)
            z = rand(Normal(mean1, 1))
            #evaluate Infected next  step - rate infection with additive noise
            S1[j, 2] = S1[j-1, 2] + (beta[2] * h  + desvest[2] * sqrt(h) * z) * (1.0 - S1[j-1, 2]) * S1[j-1, 2] - gamma * S1[j-1, 2] * h#Infected
            S1[j, 1] = 1.0 - S1[j, 2]
        end
        #saving evolution of second day (day is divided in T parts - T+1 extrems)
		S[2] = S1[T + 1,1]
		I[2] = S1[T + 1,2]
		
		#updating
		S1 = zeros(T + 1,2)
		S1[1,1] = S[2]
		S1[1,2] = I[2]
		

        bt = beta[2]
        dev = desvest[2] #fijo - cambiar

        for i in 3: Time
			for j in 2: (T + 1)
				z=rand(Normal(mean1,desv))
				#evaluate Infected next  step - rate infection with additive noise
				S1[j,2]= S1[j-1,2]+ (bt + z)*(1.0 - S1[j-1,2])*S1[j-1,2]*h - gamma*S1[j-1,2]*h;#Infected
				#evaluate Infected next  step - rate infection with proportional noise
				#S1[i+1,2]= S1[i,2] + bt*(1.0 + z)*(1.0 - S1[i,2])*S1[i,2]*h - gamma*S1[i,2]*h;#Infected
				S1[j,1]= 1.0 - S1[j,2]
			end
			#saving
			S[i] = S1[T + 1,1]
			I[i] = S1[T + 1,2]
			
			#updating values
			S1 = zeros(T + 1,3)
			S1[1,1] = S[i]
			S1[1,2] = I[i]
		    bt = beta[i]
            desvest[i]
        end
    return [S I]
end
