function calculateCO2SARF(C, N; C0=277.15)
    #constants/formula from Meinshausen et al 2020
    a1 = -2.4785e-7 #W/m2/ppm^2
    b1 = 0.00075906 #W/m2/ppm
    c1 = -0.0021492 #W/m2/ppb^0.5
    d1 = 5.2488     #W/m2
    # C0: initial CO2 conc in ppm 

    C_α_max = C0 - b1/(2*a1) # approx. = 1808 ppm
    if C > C_α_max
        α = d1 - b1^2/(4*a1)
    elseif C < C0
        α = d1
    else
        α = d1 + a1*(C-C0)^2 + b1*(C-C0)
    end

    α_N2O = c1*sqrt(N)

    return (α + α_N2O) * log(C/C0) #W/m^2
end

function calculateN2OSARF(C, N, M; N0=273.87)
    #constants/forumla from Meinshausen et al 2020
    a2 = -0.00034197 #W/m2/ppm
    b2 = 0.00025455  #W/m2/ppb
    c2 = -0.00024357 #W/m2/ppb
    d2 = 0.12173     #W/m2/ppb^0.5
    # N0: initial N2O conc 

    SARF = ( a2*sqrt(C) + b2*sqrt(N) + c2*sqrt(M) + d2 ) * ( sqrt(N) - sqrt(N0) )
    return SARF
end

function calculateCH4SARF(N, M; M0=731.41)
    #constants/forumla from Meinshausen et al 2020
    a3 = -8.9603e-5  #W/m2/ppb
    b3 = -0.00012462 #W/m2/ppb
    d3 = 0.045194    #W/m2/ppb^0.5
    # M0: initial CH4 conc

    SARF = ( a3*sqrt(M) + b3*sqrt(N) + d3 ) * ( sqrt(M) - sqrt(M0) )
    return SARF
end