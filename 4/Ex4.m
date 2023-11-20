run('Ex04_variables.m');

x0 = [0,0,0];
clockOffRe = 0;
c = s_light;
maxIterations = 10;
aed_sat = zeros(length(xyz_sat),1);

xyz_rec=x0;
for i = 1:1
    
    [azimuth, elevation, rho] = topocent(xyz_rec, xyz_sat);    
    
    xGeo = gc2gg(xyz_rec);

    tropoC = tropo_correction (xGeo(3), elevation);
    ionoC = iono_correction(xGeo(1), xGeo(2), azimuth, elevation, time_rx, ionoparams);
    
    deltaP = pr_C1 - rho - tropoC - ionoC + c*dtS;
    
    

end

