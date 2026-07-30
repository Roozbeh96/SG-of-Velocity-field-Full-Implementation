function  production_comp(obj)

    fprintf('Computing Production...\n');
    [dUdz] = gradient(mean(mean(obj.Gen_u_HRVF,3),2), obj.z(2)-obj.z(1));
    obj.production = mean(-sum((obj.Gen_u_HRVF - mean(obj.Gen_u_HRVF,3)) ...
    .*(obj.Gen_w_HRVF - mean(obj.Gen_w_HRVF,3)),3) ...
    / (size(obj.Gen_u_HRVF,3) - 1),2).*dUdz;

    obj.eta_prod = (obj.nu^3./obj.production).^0.25;

end