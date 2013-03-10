#pragma acc kernels, copy(call_result[0:nb_opt], put_result[0:nb_opt], &
#pragma acc & option_strike[0:nb_opt], stock_price[0:nb_opt], option_years[0:nb_opt]) 
{
int opt;
for(opt = 0; opt < nb_opt; opt++) {
	float sqrtT, expRT, K, d1, d2, CNDD1, CNDD2, 
		Riskfree=RISKFREE, 
		Volatility=VOLATILITY; 
	sqrtT = sqrtf(option_years[opt]);
	d1 = (logf(stock_price[opt] / option_strike[opt]) + (Riskfree +
            0.5f*Volatility*Volatility)*option_years[opt]) / (Volatility*sqrtT);
    	d2 = d1 - Volatility * sqrtT;
    
	K = 1.0f / (1.0f + 0.2316419f * fabsf(d1));
    	CNDD1 = RSQRT2PI * expf(- 0.5f*d1*d1) * (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));
	K = 1.0f / (1.0f + 0.2316419f * fabsf(d2));
    	CNDD2 = RSQRT2PI * expf(- 0.5f*d2*d2) * (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));

	expRT = expf(- Riskfree * option_years[opt]);
    	call_result[opt] = stock_price[opt] * CNDD1 - option_strike[opt] * expRT * CNDD2;
    	put_result[opt]  = option_strike[opt]*expRT*(1.0f - CNDD2) - stock_price[opt]*(1.0f - CNDD1);
	} 
}
