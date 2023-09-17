import numpy as np
import scipy.stats as si



class EuropeanOption_MC:

    def __init__(self,S,K,T,r,q,sigma,Option_Type,n,m):
        self.S = S #Spot price
        self.K = K #Strike
        self.T = T #Time to Maturity
        self.r = r #Interest rate
        self.q = q #Div
        self.sigma = sigma #vol
        self.Option_Type = Option_Type
        self.n = n
        self.m = m
        self.expected_value = 0
        self.expected_value_discounted = 0
        self.var = 0
        self.std_dev = 0




    def monte_carlo(self):
        '''
    Simple Monte Carlo option pricing algorithm for European Calls and Puts with Euler discretization
    Returns the simulated price of a European call or put with standard deviation and standard error
    associated with the simulation
            '''
        S = self.S
        K = self.K
        T =self.T  # Time to Maturity
        r = self.r   # Interest rate
        q = self.q
        sigma = self.sigma  # vol
        Option_Type = self.Option_Type
        n = self.n
        m = self.m

        @staticmethod
        def path_sim(mu, sigma, T, S, N, M):
            sims = np.zeros(M)
            dt = T / N

            for i in range(M):
                W = [0] + np.random.standard_normal(size=N)
                sims[i] = np.sum(W) * np.sqrt(dt)  # We only are concerned with the terminal value for European options

            St = S * np.exp((mu - 0.5 * sigma ** 2) * T + sigma * sims)
            return St

        paths = path_sim(r - q, sigma, T, S, n, m)


        if Option_Type == 'Call':
            paths = paths - K
            paths[paths < 0] = 0

        elif Option_Type == 'Put':
            paths = K - paths
            paths[paths < 0] = 0

        c_hat = np.mean(paths) * np.exp(-r * T)

        sd = np.sqrt(np.sum((paths * np.exp(-r * T) - c_hat) ** 2) / (m - 1))
        se = sd / np.sqrt(m)
        self.expected_value_discounted, self.var, self.std_dev = c_hat, sd, se


    def monte_carlo_anti_variates(self):
        S = self.S
        K = self.K
        T = self.T  # Time to Maturity
        r = self.r  # Interest rate
        q = self.q
        sigma = self.sigma  # vol
        Option_Type = self.Option_Type
        n = self.n
        m = self.m


        @staticmethod
        def path_sim(mu, sigma, T, S, N, M):
            sims = np.zeros(M)
            dt = T / N

            for i in range(M):
                W = [0] + np.random.standard_normal(size=N)
                sims[i] = np.sum(W) * np.sqrt(dt)  # we are concerned with only the final value of the path

            St = S * np.exp((mu - 0.5 * sigma ** 2) * T + sigma * sims)
            Sta = S * np.exp((mu - 0.5 * sigma ** 2) * T - sigma * sims)
            return np.array([St, Sta])


        paths = path_sim(r - q, sigma, T, S, n, m)

        if Option_Type == 'Call':
            paths = paths - K
            paths[paths < 0] = 0

        elif Option_Type == 'Put':
            paths = K - paths
            paths[paths < 0] = 0



        paths_mean = np.mean(paths, axis=0)
        path_mean_discount = np.mean(paths) * np.exp(-r * T)

        var = np.sqrt(np.sum((paths_mean * np.exp(-r * T) - path_mean_discount) ** 2) / (m - 1))
        std_dev = var / np.sqrt(m)
        self.expected_value_discounted, self.var, self.std_dev = path_mean_discount, var, std_dev


