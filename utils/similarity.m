
function v = similarity(mu1, spd1, mu2, spd2)
    v = exp( -( (mu1-mu2)/(spd1+spd2) )^2 );
end
function sim = possibility(mu1, spd1, mu2, spd2)
    spdmin = min(spd1, spd2);
    beta = 2 * spdmin / (spd1 + spd2);
    v = rankify(mu1, spd1, mu2, spd2);
    rtlnv = sqrt(-log(v));
    psi = beta ...
        + (1 - beta) * erf( (1/(1-beta)) *rtlnv ) ...
        - erf( rtlnv );
    sim = psi / (2 - psi);
end
