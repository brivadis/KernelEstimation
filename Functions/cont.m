function u = cont(z1, zr1, zr2, What11, What12, alpha)
    u = (1-alpha)*z1 - What11*S(zr1) - What12*S(zr2);
end