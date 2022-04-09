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

    local mp = 50
    local np = 100
    local ntotal = mp * np

    local xl = 1.0
    local yl = 2.0
    local dx = xl / mp
    local dy = yl / np

    -- 初始化多维数组
    x[1] = {}
    x[2] = {}
    vx[1] = {}
    vx[2] = {}

    for i = 1, ntotal, 1 do
        vx[1][i] = 0.0
        vx[2][i] = 0.0
        rho[i] = 1000.0
        mass[i] = dx * dy * rho[i]
        u[i] = 357.1
        itype[i] = 2 -- 淡水
        hsml[i] = dx
    end

    for i = 1, mp, 1 do
        for j = 1, np, 1 do
            local k = j + (i - 1) * np
            x[1][k] = (i - 1) * dx + dx / 2
            x[2][k] = (j - 1) * dy + dy / 2
            p[k] = rho[k] * 9.8 * (yl - x[2][k]) + 101.325e3
        end
    end

    return x, vx, mass, rho, p, u, itype, hsml, ntotal

end

-- 初始化虚粒子 @todo: 部署 Ⅱ 型虚粒子
function virt_part()
    local mp = 200
    local nvirt = 0
    local xl = 4.0
    local dx = xl / mp

    local x = {}
    local vx = {}
    local mass = {}
    local rho = {}
    local p = {}
    local u = {}
    local itype = {}
    local hsml = {}

    -- 初始化多维数组
    x[1] = {}
    x[2] = {}
    vx[1] = {}
    vx[2] = {}

    -- 上边界
    for i = 1, 2 * mp + 1 do
        nvirt = nvirt + 1
        x[1][nvirt] = (i - 1) * dx / 2
        x[2][nvirt] = xl
        vx[1][nvirt] = 0.0
        vx[2][nvirt] = 0.0
    end

    -- 下边界
    for i = 1, 2 * mp + 1 do
        for j = 1, 3, 1 do
            nvirt = nvirt + 1
            x[1][nvirt] = (i - 1) * dx / 2
            x[2][nvirt] = 0.0 - dx / 2 * (j - 1)
            vx[1][nvirt] = 0.0
            vx[2][nvirt] = 0.0
        end
    end

    -- 左边界
    for i = 1, 2 * mp - 1 do
        for j = 1, 3 do
            nvirt = nvirt + 1
            x[1][nvirt] = 0.0 - dx / 2 * (j - 1)
            x[2][nvirt] = i * dx / 2
            vx[1][nvirt] = 0.0
            vx[2][nvirt] = 0.0
        end
    end

    -- 右边界
    for i = 1, 2 * mp - 1 do
        for j = 1, 3 do
            nvirt = nvirt + 1
            x[1][nvirt] = xl + dx / 2 * (j - 1)
            x[2][nvirt] = i * dx / 2
            vx[1][nvirt] = 0.0
            vx[2][nvirt] = 0.0
        end
    end

    for i = 1, nvirt do
        rho[i] = 1000.0
        mass[i] = rho[i] * dx * dx
        p[i] = 101.325e3
        u[i] = 357.1
        itype[i] = -2 -- 虚粒子淡水
        hsml[i] = dx
    end

    return x, vx, mass, rho, p, u, itype, hsml, nvirt

end
