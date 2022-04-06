-- dofile("luanames.lua")
-- x, vx, mass, rho, p, u, itype, hsml, ntotal = input()
-- 初始化2维剪切腔体的流体实粒子
function input()
    local x = {}
    local vx = {}
    local mass = {}
    local rho = {}
    local p = {}
    local u = {}
    local itype = {}
    local hsml = {}
    local ntotal

    local m = 41
    local n = 41
    local mp = m - 1
    local np = n - 1
    local ntotal = mp * np

    local xl = 1.0e-3
    local yl = 1.0e-3
    local dx = xl / mp
    local dy = yl / np

    -- 初始化多维数组
    x[1] = {}
    x[2] = {}
    vx[1] = {}
    vx[2] = {}

    for i = 1, mp, 1 do
        for j = 1, np, 1 do
            local k = j + (i - 1) * np
            x[1][k] = (i - 1) * dx + dx / 2
            x[2][k] = (j - 1) * dy + dy / 2
        end
    end

    for i = 1, ntotal, 1 do
        vx[1][i] = 0.0
        vx[2][i] = 0.0
        rho[i] = 1000.0
        mass[i] = dx * dy * rho[i]
        p[i] = 0.0
        u[i] = 357.1
        itype[i] = 2
        hsml[i] = dx
    end

    return x, vx, mass, rho, p, u, itype, hsml, ntotal

end
