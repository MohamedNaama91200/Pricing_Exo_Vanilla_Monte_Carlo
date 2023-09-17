import numpy as np

class Barrier_Option_MC:
    def __init__(self, S, K, T, r, q, sigma, option_type, barrier_type, barrier_price, n, m):
        self.S = S  # Prix actuel de l'actif sous-jacent
        self.K = K  # Prix d'exercice
        self.T = T  # Temps à l'échéance
        self.r = r  # Taux d'intérêt
        self.q = q  # Taux de dividende
        self.sigma = sigma  # Volatilité
        self.option_type = option_type  # Type d'option (Call ou Put)
        self.barrier_type = barrier_type  # Type de barrière (Up-and-In, Up-and-Out, Down-and-In, Down-and-Out)
        self.barrier_price = barrier_price  # Prix de la barrière
        self.n = n  # Nombre de pas dans le modèle de Monte Carlo
        self.m = m  # Nombre de simulations
        self.expected_value_discounted = 0
        self.var = 0
        self.std_dev = 0

    def path_simulation(self, mu, sigma, T, S, N, M):
        sims = np.zeros((M, N))
        dt = T / N

        for i in range(M):
            W = np.random.standard_normal(size=N)
            sims[i, 0] = S
            for j in range(1, N):
                sims[i, j] = sims[i, j - 1] * np.exp((mu - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * W[j])

        return sims

    def monte_carlo(self):
        paths = self.path_simulation(self.r - self.q, self.sigma, self.T, self.S, self.n, self.m)


        if self.barrier_type == 'Up-and-In':
            barrier_condition = np.all(paths < self.barrier_price, axis=1)
        elif self.barrier_type == 'Up-and-Out':
            barrier_condition = np.any(paths > self.barrier_price, axis=1)
        elif self.barrier_type == 'Down-and-In':
            barrier_condition = np.all(paths > self.barrier_price, axis=1)
        elif self.barrier_type == 'Down-and-Out':
            barrier_condition = np.any(paths < self.barrier_price, axis=1)



        if self.option_type == 'Call':
            if barrier_condition.any():
                payoff = np.maximum(paths[:, -1] - self.K, 0)
                payoff[barrier_condition] = 0
            else:
                payoff = np.maximum(paths[:, -1] - self.K, 0)
        elif self.option_type == 'Put':
            if barrier_condition.any():
                payoff = np.maximum(self.K - paths[:, -1], 0)
                payoff[barrier_condition] = 0

            else:
                payoff = np.maximum(self.K - paths[:, -1], 0)


        c_hat = np.mean(payoff) * np.exp(-self.r * self.T)
        sd = np.std(payoff) / np.sqrt(self.m)
        self.expected_value_discounted, self.var, self.std_dev = c_hat, sd, sd

S, K, T, r, q, sigma, Option_Type,Barrier_Type,Barriere_Price, n, m = 100,100,1,0.05,0,0.3,"Put","Down-and-In",80, 5,10
Barriere = Barrier_Option_MC(S, K, T, r, q, sigma, Option_Type,Barrier_Type,Barriere_Price, n, m)
Barriere.monte_carlo()
