#  Function SIS let us obtain  values for:
#  S_n,I_n
#  T is the steps number.

using Random, Distributions
include("./create_forward_beta_estimator.jl")
include("./create_forward_sigma_estimator.jl")

function SDE_SIS(T_extrem::Int64, divitions::Int64, Num_subdivitions::Int64, mean1::Float64, desvest::Any, gamma::Float64, beta::Any, I0::Float64)
    #SIS(T,Time,h,mean,desvest,gamma,beta,I0)
    Time_interval = range(0.00, T_extrem, divitions * Num_subdivitions) # Num_subdivitions
    T = Num_subdivitions
    tj = (Time_interval[2] - Time_interval[1])
    h = tj
    Time =  divitions
    #initial disease rate
    
    S1 = zeros(T + 1, 2)#auxiliar subdivitions steps

    #create_vectors of returns
    S = zeros(Time) # Suceptibles in divitions steps 
    I = zeros(Time) # Infected in divitions steps
    B = zeros(Time)
    Sigma = zeros(Time)
    MLE_beta = zeros(Time) # Maximum likelihood beta estimator 

    #Initialization of I_0 S_0 
    S1[1, 2] = I0[1] #I_0
    S1[1, 1] = 1.0 - I0[1] #S_0
    S[1] = S1[1, 1]
    I[1] = S1[1, 2]
    # SDE SIS FIRST EVOLUTION
    sum_MLE_beta = 0.0
    sum_numerator_sigma_estimator = 0.0
    sum_denominator_sigma_estimator = 0.0
    for j in 2: (T+1)
        z = rand(Normal(mean1, 1)) 
        #evaluate Infected next  step j - rate infection with additive noise evaluating with step j-1 data
        S1[j, 2] = S1[j-1, 2] + (beta[j-1] * h + desvest[j - 1]*sqrt(h)*z) * (1.0 - S1[j-1, 2]) * S1[j-1, 2]  - gamma * S1[j-1, 2] * h#Infected
        S1[j, 1] = 1.0 - S1[j, 2] #Suceptibles
        #S1[j,3] = create_beta_estimator(S1[:,2],gamma,tj,j,j)
        #S1[j,4] = create_sigma_estimator(S1[:,2],j)
        sum_MLE_beta = sum_MLE_beta + (S1[j, 2] - S1[j-1, 2]) / ( S1[j - 1, 2] * (1 - S1[j-1, 2] ) ) + gamma * h / ( (1 - S1[j - 1, 2] ) )       
        sum_numerator_sigma_estimator = sum_numerator_sigma_estimator + (S1[j, 2] - S1[j-1, 2])^2 #( S1[j, 2] * (1 - S1[j, 2] ) - S1[j - 1, 2] * (1 - S1[j-1, 2] ) )^2 # 
        sum_denominator_sigma_estimator = sum_denominator_sigma_estimator + ( S1[j - 1, 2] * (1 - S1[j-1, 2] ) )^2 * h # + ( S1[j, 2] * (1 - S1[j, 2] ) )^2
    end
    #saving the final values
    S[2] = S1[T + 1, 1]
    I[2] = S1[T + 1, 2]
    #B[2] = mean(S1[:,3])
    Sigma[1] = sqrt( sum_numerator_sigma_estimator * (sum_denominator_sigma_estimator )^(-1)) #* (Time_interval[T]-Time_interval[1])
    MLE_beta[1] = sum_MLE_beta

        S1 = zeros(T + 1,2)
        S1[1,1] = S[2] # redefine initial condition
        S1[1,2] = 1 - S1[1,1]
        sum_MLE_beta1 = 0.0
        sum_numerator_sigma_estimator1 = 0.0
        sum_denominator_sigma_estimator1 = 0.0
        for j in 2:(T+1)
            z = rand(Normal(mean1, 1))
            #evaluate Infected next  step - rate infection with additive noise
            S1[j, 2] = S1[j-1, 2] + (beta[T + j-1] * h + desvest[T + (j-1) ]*sqrt(h)*z) * (1.0 - S1[j-1, 2]) * S1[j-1, 2] - gamma * S1[j-1, 2] * h#Infected
            S1[j, 1] = 1.0 - S1[j, 2]
            #S1[j,3] = create_beta_estimator(S1[:,2],gamma,tj,j,j)
            #S1[j,4] = create_sigma_estimator(S1[:,2],j)
            sum_MLE_beta1 = sum_MLE_beta1 + (S1[j, 2] - S1[j-1, 2]) / ( S1[j - 1, 2] * (1 - S1[j-1, 2] ) ) + gamma * h / ( (1 - S1[j - 1, 2] ) ) 
            sum_numerator_sigma_estimator1 = sum_numerator_sigma_estimator1 + (S1[j, 2] - S1[j-1, 2])^2  
            sum_denominator_sigma_estimator1 = sum_denominator_sigma_estimator1 + ( S1[j - 1, 2] * (1 - S1[j-1, 2] ) )^2 * h 
        end
        #saving evolution of second day (day is divided in T parts - T+1 extrems)
		S[3] = S1[T + 1,1]
		I[3] = S1[T + 1,2]
        #B[3] = mean(S1[:,3])
        
        MLE_beta[2] = MLE_beta[1] + sum_MLE_beta1
        MLE_beta[1] = MLE_beta[1]*(Time_interval[T]-Time_interval[1])^(-1)
        Sigma[2] = sqrt( sum_numerator_sigma_estimator1 * (sum_denominator_sigma_estimator1 )^(-1))


		#updating
		S1 = zeros(T + 1,2)
		S1[1,1] = S[3]
		S1[1,2] = I[3]
		
        
        for i in 3: Time
            sum_MLE_betai = 0.0
            sum_numerator_sigma_estimatori = 0.0
            sum_denominator_sigma_estimatori = 0.0
			for j in 2: (T + 1)
				z=rand(Normal(mean1,1))
				#evaluate Infected next  step - rate infection with additive noise
				S1[j,2]= S1[j-1,2]+ (beta[(i-1)*T + j-1] * h + desvest[(i-1)*T + j-1]*sqrt(h)*z)*(1.0 - S1[j-1,2])*S1[j-1,2] - gamma*S1[j-1,2]*h;#Infected
				S1[j,1]= 1.0 - S1[j,2]
                #S1[j,3] = create_beta_estimator(S1[:,2],gamma,tj,j,j)
                #S1[j,4] = create_sigma_estimator(S1[:,2],j)
                sum_MLE_betai = sum_MLE_betai + (S1[j, 2] - S1[j-1, 2]) / ( S1[j - 1, 2] * (1 - S1[j-1, 2] ) ) + gamma * h / ( (1 - S1[j - 1, 2] ) ) 
                sum_numerator_sigma_estimatori = sum_numerator_sigma_estimatori + (S1[j, 2] - S1[j-1, 2])^2  
            sum_denominator_sigma_estimatori = sum_denominator_sigma_estimatori + ( S1[j - 1, 2] * (1 - S1[j-1, 2] ) )^2 * h 
            end
			#saving
			S[i] = S1[T + 1,1]
			I[i] = S1[T + 1,2]
            #B[i] = mean(S1[:,3])
            
			MLE_beta[i] = MLE_beta[i - 1] + sum_MLE_betai
            sum_MLE_betai = 0.0
            MLE_beta[i - 1] = MLE_beta[i - 1] * (Time_interval[(i-1)*T] - Time_interval[1])^(-1)
            Sigma[i] = sqrt( sum_numerator_sigma_estimatori * (sum_denominator_sigma_estimatori )^(-1))
            #updating values
			S1 = zeros(T + 1,2)
			S1[1,1] = S[i]
			S1[1,2] = I[i]
        end 
    
    # Get MLE_beta_estimator
    #for i in 3: Time

    #end

    MLE_beta[Time] = MLE_beta[Time] * (T_extrem)^(-1)  
    return [S I MLE_beta Sigma ]
end
