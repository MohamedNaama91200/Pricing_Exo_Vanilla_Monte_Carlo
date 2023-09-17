import numpy as np

class Asian_Option_MC:
    def __init__(self, S, K, T, r, q, sigma, option_type, n, m):
        self.S = S  # Prix actuel de l'actif sous-jacent
        self.K = K  # Prix d'exercice
        self.T = T  # Temps à l'échéance
        self.r = r  # Taux d'intérêt
        self.q = q  # Taux de dividende
        self.sigma = sigma  # Volatilité
        self.option_type = option_type  # Type d'option (Call ou Put)
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
        average_paths = np.mean(paths, axis=1)

        if self.option_type == 'Call':
            payoff = average_paths - self.K
            payoff[payoff < 0] = 0
        elif self.option_type == 'Put':
            payoff = self.K - average_paths
            payoff[payoff < 0] = 0

        c_hat = np.mean(payoff) * np.exp(-self.r * self.T)
        sd = np.std(payoff) / np.sqrt(self.m)
        self.expected_value_discounted, self.var, self.std_dev = c_hat, sd, sd

    def monte_carlo_anti_variates(self):
        paths = self.path_simulation(self.r - self.q, self.sigma, self.T, self.S, self.n, self.m)
        average_paths = np.mean(paths, axis=1)

        if self.option_type == 'Call':
            payoff = average_paths - self.K
            payoff[payoff < 0] = 0
        elif self.option_type == 'Put':
            payoff = self.K - average_paths
            payoff[payoff < 0] = 0

        paths_mean = np.mean(payoff)
        path_mean_discount = paths_mean * np.exp(-self.r * self.T)

        var = np.sqrt(np.sum((paths_mean * np.exp(-self.r * self.T) - path_mean_discount) ** 2) / (self.m - 1))
        std_dev = var / np.sqrt(self.m)
        self.expected_value_discounted, self.var, self.std_dev = path_mean_discount, var, std_dev
