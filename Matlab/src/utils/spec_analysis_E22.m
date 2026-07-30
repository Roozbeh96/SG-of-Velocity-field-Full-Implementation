function spec_analysis_E22(obj)
    %{
    More info about how to compute the structure function in:
    https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/
    article/stochastic-modal-velocity-field-in-roughwall-turbulence/
    800D97E60609F4FE4D92847CCDB1C7A0
    %}
    fprintf('Spectral Analysis E22...\n');
    wprime = obj.Gen_w_HRVF - mean(obj.Gen_w_HRVF,3);
    concat_wprime = reshape(wprime,size(wprime, 1), []);
    ksmax=floor(2*pi/(obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1)));
    
    win15del=floor(15*obj.delta/(obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1)));
    window=hann(win15del);
    NFFT=2^(nextpow2(win15del));
    [PSD_k,k]=pwelch(concat_wprime',window,NFFT/2,win15del,ksmax); %PSD_f[m^2/s]
    Apk = trapz(k,PSD_k);
    Q = mean(var(obj.Gen_w_HRVF,0,3),2)./Apk';  %%[m^3/s^2]=[m^2/s]*[m^2/s^2]*1/[m/s]
    obj.E22=(PSD_k.*Q')';
end