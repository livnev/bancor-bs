import numpy as np

#ETH_SUPPLY = 100e6
DAI_SUPPLY = 50e6
BNT_SUPPLY = 75e6

class Reserve:
    def __init__(self, bal_dai=270000, bal_bnt=55000, initial_daibnt=92000, fee=0.001, connector_weight=500000):
        self.bal_dai = bal_dai
        self.bal_bnt = bal_bnt
        self.connector_weight = connector_weight
        self.fee = fee
        self.daibnt_supply = initial_daibnt

    def balance_check(self):
        assert self.bal_dai >= 0
        assert self.bal_bnt >= 0
        
    def buy(self, wad):
        # buy BNT for DAI, going through DAIBNT token
        # wad : dai

        # first buys DAIBNT for DAI
        ddai, ddaibnt = self.buy_daibnt_dai(wad)
        # then sells DAIBNT for BNT
        dbnt, _ = self.sell_daibnt_bnt(ddaibnt)
        
        self.bal_dai -= ddai
        self.bal_bnt -= dbnt
        self.balance_check()
        return (ddai, dbnt)

    def sell(self, wad):
        # sell BNT for DAI, going through DAIBNT token
        # wad : BNT

        # first buys DAIBNT for BNT
        dbnt, ddaibnt = self.buy_daibnt_bnt(wad)
        # then sells DAIBNT for DAI
        ddai, _ = self.sell_daibnt_dai(ddaibnt)
        
        self.bal_dai -= ddai
        self.bal_bnt -= dbnt
        self.balance_check()
        return (ddai, dbnt)

    def buy_daibnt_dai(self, wad):
        buy_daibnt = self.get_purchase_return_dai(wad)
        self.daibnt_supply += buy_daibnt
        return (-wad, buy_daibnt)
        
    def buy_daibnt_bnt(self, wad):
        buy_daibnt = self.get_purchase_return_bnt(wad)
        self.daibnt_supply += buy_daibnt
        return (-wad, buy_daibnt)
        
    def sell_daibnt_dai(self, wad):
        buy_dai = self.get_sale_return_dai(wad)
        self.daibnt_supply -= wad
        return (buy_dai, -wad)
    
    def sell_daibnt_bnt(self, wad):
        buy_bnt = self.get_sale_return_bnt(wad)
        self.daibnt_supply -= wad
        return (buy_bnt, -wad)
        
    def get_mid_price(self):
        # hack
        return 0.5 * (self.get_buy_price(1e-5) + self.get_sell_price(1e-5))
        
    def get_buy_price(self, wad):
        # gets BNT:DAI buy price
        # save initial balances
        _bal_dai, _bal_bnt, _daibnt_supply = self.bal_dai, self.bal_bnt, self.daibnt_supply
        _, buy_daibnt = self.buy_daibnt_dai(wad)        
        buy_bnt, _ = self.sell_daibnt_bnt(buy_daibnt)
        price = wad / buy_bnt
        # revert balances
        self.bal_dai, self.bal_bnt, self.daibnt_supply = _bal_dai, _bal_bnt, _daibnt_supply
        return price
        
    def get_sell_price(self, wad):
        # gets BNT:DAI sell price
        # save initial balances
        _bal_dai, _bal_bnt, _daibnt_supply = self.bal_dai, self.bal_bnt, self.daibnt_supply
        _, buy_daibnt = self.buy_daibnt_bnt(wad)
        buy_dai, _ = self.sell_daibnt_dai(buy_daibnt)
        price = buy_dai / wad
        # revert balances
        self.bal_dai, self.bal_bnt, self.daibnt_supply = _bal_dai, _bal_bnt, _daibnt_supply
        return price

    # these formulas are taken directly from the Bancor contract at:
    # https://github.com/bancorprotocol/contracts/blob/master/solidity/contracts/BancorFormula.sol
    def get_purchase_return_dai(self, wad):
        # wad : dai
        return self.daibnt_supply * ((1 + wad / self.bal_dai) ** (self.connector_weight / 1000000) - 1) * (1 - self.fee)

    def get_purchase_return_bnt(self, wad):
        # wad : bnt
        return self.daibnt_supply * ((1 + wad / self.bal_bnt) ** (self.connector_weight / 1000000) - 1) * (1 - self.fee)

    def get_sale_return_dai(self, wad):
        # wad : daibnt
        return self.bal_dai * (1 - (1 - wad / self.daibnt_supply) ** (1 / (self.connector_weight / 1000000))) * (1 - self.fee)
    
    def get_sale_return_bnt(self, wad):
        # wad : daibnt
        return self.bal_bnt * (1 - (1 - wad / self.daibnt_supply) ** (1 / (self.connector_weight / 1000000))) * (1 - self.fee)

    # trick to calculate the maximum amount that can be bought or sold at a certain limit price
    def get_buy_amount(self, limit, num_tries=10):
        try_wad = self.bal_dai/16
        for i in range(num_tries):
            if self.get_buy_price(try_wad) < limit:
                return try_wad
            try_wad /= 2
        return 0.

    def get_sell_amount(self, limit, num_tries=10):
        try_wad = self.bal_bnt/16
        for i in range(100):
            if self.get_sell_price(try_wad) > limit:
                return try_wad
            try_wad /= 2
        return 0.


def sim_arb_trading(reserve, prices):
    # prices is an array of prices
    # at each time step the trader will trade with the reserve
    # such that it buys for less than prices[i]
    # and sells for more than prices[i]
    # in particular, it will try to trade at the midpoint between the market price
    # and the price offered by the Bancor reserve

    # can also start with a positive dai balance
    trader_dai = 0
    trader_bnt = 0

    bals_dai = []
    bals_bnt = []
    
    for price in prices:
        bals_dai.append(reserve.bal_dai)
        bals_bnt.append(reserve.bal_bnt)
        mid_price = reserve.get_mid_price()
        if mid_price > price:
            wad = reserve.get_sell_amount(0.5*(price+mid_price))
            if wad == 0:
                continue
            # get the bnt from the market:
            trader_dai -= wad * price
            trader_bnt += wad
            # trade with the reserve
            ddai, dbnt = reserve.sell(wad)
            trader_dai += ddai
            trader_bnt += dbnt
        elif mid_price < price:
            wad = reserve.get_buy_amount(0.5*(price+mid_price))
            if wad == 0:
                continue
            # trade with the reserve
            ddai, dbnt = reserve.buy(wad)
            trader_dai += ddai
            trader_bnt += dbnt
            # now dump the bnt
            trader_dai += trader_bnt * price
            trader_bnt -= trader_bnt
        assert trader_bnt == 0
        assert trader_dai >= 0

    return trader_dai, trader_bnt, bals_dai, bals_bnt

def simulate_geometric_wiener_old(T, start_price, mu=0., sigma=.01):
    # e.g. can model price martingale by a
    # geometric process with normal steps
    # i.e. brownian motion
    steps  = 1 + np.random.normal(size=T, loc=mu, scale=sigma)
    return start_price * np.cumprod(steps)

def simulate_geometric_wiener(T, start_price, mu=0., sigma=.01, dt=0.1):
    steps = np.exp((mu - sigma**2/2)*dt) * np.exp(sigma * np.random.normal(0, np.sqrt(dt),(1, T)))
    return start_price * np.cumprod(steps)

def simulate_garch_prices(T, start_price, a0=0.2, a1=0.5, b1=0.3, mu=0., sigma=.01):
    # simulates a GARCH(1, 1) process
    w = np.random.normal(size=T, loc=mu, scale=sigma)
    eps = np.zeros_like(w)
    sigsq = np.zeros_like(w)
    for i in range(1, T):
        sigsq[i] = a0 + a1*(eps[i-1]**2) + b1*sigsq[i-1]
        eps[i] = w[i] * np.sqrt(sigsq[i])
    return start_price * np.cumprod(1 + eps)
