function y = peso(n, snr, var)
  
nueva_desv = sqrt(var);  
fdp_ruido = normpdf(n, 0, sqrt(1/(2*snr)));
fdp_ruido_estrella = normpdf(n, 0, nueva_desv);

y = (fdp_ruido./fdp_ruido_estrella)';
 
 