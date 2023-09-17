import numpy
import scipy.stats as si
from European_Option_MC import *
from Asian_Option_MC import *
from Barrier_Option_MC import *

def main() :

    S,K,T,r,q,sigma,Option_Type,n,m = get_user_input()
    MC = EuropeanOption_MC(S, K, T, r, q, sigma, Option_Type, n, m)
    MC.monte_carlo()
    price = MC.expected_value_discounted
    std_dev = MC.std_dev
    print(f"Prix de l'option par MC : {price}")
    print(f"Ecart-Type :{std_dev}")

    MC_Anti_Variates = EuropeanOption_MC(S,K,T,r,q,sigma,Option_Type,n,m)
    MC_Anti_Variates.monte_carlo_anti_variates()
    price = MC_Anti_Variates.expected_value_discounted
    std_dev = MC_Anti_Variates.std_dev
    print(f"Prix de l'option par MC Anti-Variates : {price}")
    print(f"Ecart-Type :{std_dev}")

    S, K, T, r, q, sigma, Option_Type,Barrier_Type, n, m = get_user_input()
    MC = Asian_Option_MC(S, K, T, r, q, sigma, Option_Type, n, m)
    MC.monte_carlo()
    MC.monte_carlo_anti_variates()
    price = MC.expected_value_discounted
    std_dev = MC.std_dev
    print(f"Prix de l'option asiatique  par MC : {price}")
    print(f"Ecart-Type :{std_dev}")
    print(f"Prix de l'option asiatique  par MC AV : {price}")
    print(f"Ecart-Type :{std_dev}")

    MC_Anti_Variates = EuropeanOption_MC(S, K, T, r, q, sigma, Option_Type, n, m)
    MC_Anti_Variates.monte_carlo_anti_variates()
    price = MC_Anti_Variates.expected_value_discounted
    std_dev = MC_Anti_Variates.std_dev
    price_av = MC_Anti_Variates.expected_value_discounted
    std_dev_av = MC_Anti_Variates.std_dev
    print(f"Prix de l'option asiatique par MC Anti-Variates : {price_av}")
    print(f"Ecart-Type :{std_dev_av}")

    S, K, T, r, q, sigma, Option_Type, Barrier_Type,Barrier_Price ,n, m = get_user_input()
    MC = Barrier_Option_MC(S, K, T, r, q, sigma, Option_Type, Barrier_Type,Barrier_Price,n, m)
    MC.monte_carlo()
    price = MC.expected_value_discounted
    std_dev = MC.std_dev
    print(f"Prix de l'option  par MC : {price}")
    print(f"Ecart-Type :{std_dev}")





def get_user_input():
    S = float(input("Prix actuel de l'actif sous-jacent (S) : "))
    K = float(input("Prix d'exercice de l'option (K) : "))
    T = float(input("Temps à l'échéance (T) : "))
    r = float(input("Taux d'intérêt (r) : "))
    q = float(input("Taux de dividende (q) : "))
    sigma = float(input("Volatilité (sigma) : "))
    Option_Type = input("Type d'option (Call ou Put) : ")
    Barrier_Type = input("Type d'option barrières (Up-and-In,Up-and-Out,Down-and-Out,Down-and-In) : ")
    Barrier_Price = float(input("Prix de la barrière :"))

    n = int(input("Nombre de pas dans le modèle de Monte Carlo (n) : "))
    m = int(input("Nombre de simulations (m) : "))

    return S, K, T, r, q, sigma, Option_Type,Barrier_Type,Barrier_Price, n, m

if __name__ == "__main__":
    main()







