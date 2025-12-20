#include "itensor/all.h"
#include <random>
#include <complex>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>   // 用于保存 counts
#include <cstddef>   // std::size_t

using namespace itensor;

// ====================== Haar 随机 U(4) 生成 =======================

void
haarRandomU4(std::mt19937 &rng,
             std::complex<double> U[4][4])
{
    using C = std::complex<double>;
    std::normal_distribution<double> ndist(0.0,1.0);

    // 4x4 复高斯矩阵
    C Z[4][4];
    for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
    {
        double x = ndist(rng);
        double y = ndist(rng);
        Z[i][j] = C(x,y);
    }

    // 列 Gram-Schmidt -> Haar 随机 U(4)
    C Q[4][4];
    for(int j = 0; j < 4; ++j)
    {
        C v[4];
        for(int i = 0; i < 4; ++i) v[i] = Z[i][j];

        for(int k = 0; k < j; ++k)
        {
            C proj = 0.0;
            for(int i = 0; i < 4; ++i)
            {
                proj += std::conj(Q[i][k]) * v[i];
            }
            for(int i = 0; i < 4; ++i)
            {
                v[i] -= proj * Q[i][k];
            }
        }

        double norm2 = 0.0;
        for(int i = 0; i < 4; ++i) norm2 += std::norm(v[i]);
        double norm = std::sqrt(norm2);

        for(int i = 0; i < 4; ++i)
        {
            Q[i][j] = v[i] / norm;
        }
    }

    for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
    {
        U[i][j] = Q[i][j];
    }
}

// ====================== MPS 上作用两比特门 =======================

// 作用 gate G 在 (i,i+1) 上
// G 索引结构: (s1,s2,s1',s2')
void
applyTwoSiteGateMPS(MPS &psi,
                    ITensor const& G,
                    int i,
                    int maxDim)
{
    psi.position(i);

    auto wf = psi(i) * psi(i+1);

    // 更推荐：明确门作用方向，避免隐式匹配
    // 让门把 (s_i,s_{i+1}) -> (s_i',s_{i+1}')
    wf = G * wf;

    wf.noPrime();

    auto s1 = siteIndex(psi,i);
    auto li = leftLinkIndex(psi,i); // i==1 时可能不存在

    IndexSet left = (li ? IndexSet(li,s1) : IndexSet(s1));

    auto [U,S,V] = svd(wf, left,
                       Args("Cutoff=",1E-12,"MaxDim=",maxDim));

    psi.setA(i,U);
    psi.setA(i+1,S*V);

}

// 构造 Haar 随机两比特门并作用在 (i,i+1)
void
applyRandomTwoQubitGateMPS(MPS &psi,
                           int i,
                           std::mt19937 &rng,
                           int maxDim)
{
    using C = std::complex<double>;
    C U4[4][4];
    haarRandomU4(rng,U4);

    auto s1  = siteIndex(psi,i);
    auto s2  = siteIndex(psi,i+1);
    auto s1p = prime(s1);
    auto s2p = prime(s2);

    ITensor G(s1,s2,s1p,s2p);

    for(int in1=0; in1<2; ++in1)
    for(int in2=0; in2<2; ++in2)
    for(int out1=0; out1<2; ++out1)
    for(int out2=0; out2<2; ++out2)
    {
        int col = in1*2 + in2;
        int row = out1*2 + out2;

        auto inIdx1  = (in1  == 0 ? s1(1)  : s1(2));
        auto inIdx2  = (in2  == 0 ? s2(1)  : s2(2));
        auto outIdx1 = (out1 == 0 ? s1p(1) : s1p(2));
        auto outIdx2 = (out2 == 0 ? s2p(1) : s2p(2));

        G.set(inIdx1,inIdx2,outIdx1,outIdx2,U4[row][col]);
    }

    applyTwoSiteGateMPS(psi,G,i,maxDim);
}

// ====================== TFIM 基态 (DMRG) =======================
// H = -J * sum sigma^z_i sigma^z_{i+1} - h * sum sigma^x_i
// 注意：SpinHalf 的 Sz/Sx 是自旋算符 (1/2)*sigma，所以要乘因子：
// sigma^z = 2*Sz, sigma^x = 2*Sx
// => sigma^z sigma^z = 4*Sz*Sz, sigma^x = 2*Sx
MPS
makeTFIMGroundState(const SiteSet &sites, double J, double h)
{
    int N = length(sites);

    AutoMPO ampo(sites);

    for(int i = 1; i <= N-1; ++i)
    {
        ampo += -4.0*J,"Sz",i,"Sz",i+1;   // -J * sigma^z sigma^z
    }
    for(int i = 1; i <= N; ++i)
    {
        ampo += -2.0*h,"Sx",i;           // -h * sigma^x
    }

    auto H = toMPO(ampo);

    // 初始 MPS（随便给个简单直积态）
    InitState init(sites);
    for(int i = 1; i <= N; ++i) init.set(i,"Dn");
    MPS psi0(init);

    Sweeps sweeps(8);
    sweeps.maxdim() = 20, 40, 80, 120, 200, 200, 200, 200;
    sweeps.cutoff() = 1E-12;

    auto energy = dmrg(H, psi0, sweeps, Args("Quiet=",true));
    std::ostringstream oss;
    oss << "MPS/ising_N" << N
        << "_J" << std::fixed << std::setprecision(2) << J
        << "_h" << std::fixed << std::setprecision(2) << h
        << ".mps";

    std::string fname = oss.str();
    writeToFile(fname, psi0);

    // dmrg 会把 psi0 就地更新成基态
    return psi0;
}


// ====================== GHZ 态构造（MPS） =======================

void
applyHadamardMPS(MPS &psi,
                 const SiteSet &sites,
                 int i)
{
    psi.position(i);
    auto s  = sites(i);
    auto sp = prime(s);

    ITensor H(sp,s);

    auto dn  = s(1);
    auto up  = s(2);
    auto dnp = sp(1);
    auto upp = sp(2);

    double isrt2 = 1.0/std::sqrt(2.0);

    H.set(dnp,dn,isrt2);
    H.set(dnp,up,isrt2);
    H.set(upp,dn,isrt2);
    H.set(upp,up,-isrt2);

    auto Ai = psi.A(i);
    Ai = H * Ai;
    Ai.noPrime();
    psi.setA(i,Ai);
}

// CNOT 控制 i -> 受控 i+1
void
applyCNOTMPS(MPS &psi,
             const SiteSet &sites,
             int i,
             int maxDim)
{
    auto j = i+1;
    auto s1  = sites(i);
    auto s2  = sites(j);
    auto s1p = prime(s1);
    auto s2p = prime(s2);

    ITensor G(s1,s2,s1p,s2p);

    auto dn1  = s1(1);
    auto up1  = s1(2);
    auto dn2  = s2(1);
    auto up2  = s2(2);

    auto dn1p = s1p(1);
    auto up1p = s1p(2);
    auto dn2p = s2p(1);
    auto up2p = s2p(2);

    // |00> -> |00>
    G.set(dn1,dn2,dn1p,dn2p,1.0);
    // |01> -> |01>
    G.set(dn1,up2,dn1p,up2p,1.0);
    // |10> -> |11>
    G.set(up1,dn2,up1p,up2p,1.0);
    // |11> -> |10>
    G.set(up1,up2,up1p,dn2p,1.0);

    applyTwoSiteGateMPS(psi,G,i,maxDim);
}

// |00..0> -> GHZ: (|00..0> + |11..1>)/sqrt(2)
void
makeGHZMPS(MPS &psi,
           const SiteSet &sites,
           int maxDim)
{
    int N = length(sites);

    applyHadamardMPS(psi,sites,1);
    for(int i = 1; i < N; ++i)
    {
        applyCNOTMPS(psi,sites,i,maxDim);
    }
}

// ====================== brickwork 浅层线路（MPS） =======================

void
applyBrickworkCircuit(MPS &psi,
                      int N,
                      int depth,
                      std::mt19937 &rng,
                      int maxDim)
{

    for(int layer = 0; layer < depth; ++layer)
    {
        if(layer % 2 == 0)
        {
            // 偶数层：作用 (1,2), (3,4), ...
            for(int i = 1; i <= N-1; i += 2)
            {
                applyRandomTwoQubitGateMPS(psi,i,rng,maxDim);
            }
        }
        else
        {
            // 奇数层：作用 (2,3), (4,5), ...
            for(int i = 2; i <= N-1; i += 2)
            {
                applyRandomTwoQubitGateMPS(psi,i,rng,maxDim);
            }
        }
    }
}

// ====================== 计算完整概率分布 =======================

// 把 MPS 收缩成一个 wavefunction ITensor（只带所有 site index）
ITensor
makeFullWavefunction(const MPS &psi)
{
    int N = length(psi);

    auto wf = psi(1);
    for(int i = 2; i <= N; ++i)
    {
        wf *= psi(i);
    }
    wf.noPrime();
    return wf;
}


// std::vector<double> computeProbDistribution(const MPS &psi, const SiteSet &sites) {
//     int N = length(sites); 
//     auto wf = makeFullWavefunction(psi); 
//     std::size_t dim = (std::size_t(1) << N); 
//     std::vector<double> probs(dim,0.0); 
//     for(std::size_t x = 0; x < dim; ++x) { 
//         ITensor ampT = wf; 
//         // 把 wf 投影到 |x> = |b1 b2 ... bN> 
//         for(int i = 1; i <= N; ++i) { 
//             int bit = ( (x >> (i-1)) & 1 ); 
//             auto s = sites(i); ITensor ket(s); 
//             ket.set( bit==0 ? s(1) : s(2), 1.0 ); ampT *= ket; 
//         } 
//         std::complex<double> amp = eltC(ampT); 
//         probs[x] = std::norm(amp); 
//     } 
//     // 归一化 
//     double Z = 0.0; 
//     for(auto p : probs) Z += p; 
//     if(Z > 0.0) { 
//         for(auto &p : probs) p /= Z; 
//     } 
//     return probs; 
// }


std::vector<double>
computeProbDistribution(const MPS &psi)
{
    const int N = length(psi);
    const std::size_t dim = std::size_t(1) << N;

    // 复制一份 MPS，方便调整规范
    auto psiC = psi;
    psiC.position(1); // 左规范化，数值和性能都更好

    std::vector<double> probs(dim, 0.0);

    std::vector<ITensor> env_cur, env_next;
    env_cur.reserve(1);
    env_cur.emplace_back(1.0); // 空前缀 = 标量 1

    for(int i = 1; i <= N; ++i)
    {
        env_next.clear();
        env_next.reserve(env_cur.size() * 2);

        // 关键：从 psiC 取 site index（而不是外部 sites）
        auto s  = siteIndex(psiC, i);
        auto Ai = psiC.A(i);  // site tensor

        // bra |0> 和 |1>
        ITensor bra0(s), bra1(s);
        bra0.set(s(1), 1.0);
        bra1.set(s(2), 1.0);

        for(const auto &L : env_cur)
        {
            // prefix + 0
            env_next.emplace_back(L * Ai * bra0);

            // prefix + 1
            env_next.emplace_back(L * Ai * bra1);
        }

        env_cur.swap(env_next);
    }

    for(std::size_t x = 0; x < dim; ++x)
    {
        std::complex<double> amp = eltC(env_cur[x]);
        probs[x] = std::norm(amp);
    }

    double Z = 0.0;
    for(double p : probs) Z += p;
    if(Z > 0.0)
    {
        for(auto &p : probs) p /= Z;
    }

    return probs;
}



// ====================== depolarization on probabilities =======================
// p_noisy(x) = (1-p)*p(x) + p/2^N
void applyDepolarizationProb(std::vector<double> &probs, double p)
{
    std::size_t dim = probs.size();
    double uni = 1.0 / double(dim);

    // 直接在 probs 上进行修改
    for (std::size_t s = 0; s < dim; ++s)
    {
        probs[s] = (1.0 - p) * probs[s] + p * uni;
    }

    double Z = 0.0;
    for (std::size_t s = 0; s < dim; ++s) Z += probs[s];
    for (std::size_t s = 0; s < dim; ++s) probs[s] /= Z;
}



unsigned long long choose3(std::size_t M)
{
    if(M < 3) return 0ULL;
    unsigned long long m = static_cast<unsigned long long>(M);
    return m*(m-1)*(m-2)/6;
}

unsigned long long choose2(std::size_t M)
{
    if (M < 2) return 0ULL;
    unsigned long long m = static_cast<unsigned long long>(M);
    return m * (m - 1) / 2;
}

// ====================== main =======================

int
main(int argc, char* argv[])
{
    // 默认参数
    int N        = 10;   // qubit 数
    int depth    = 4;    // brickwork 深度
    int nSamples = 50;   // 采样次数
    int maxDim   = 300;  // 最大 bond dimension
    double dprob = 0.0;
    double h = 1.0;
    int experiment_id = 0;  // 默认 0

    if(argc > 1) N        = std::stoi(argv[1]);
    if(argc > 2) depth    = std::stoi(argv[2]);
    if(argc > 3) nSamples = std::stoi(argv[3]);
    if(argc > 4) maxDim   = std::stoi(argv[4]);
    if(argc > 5) dprob    = std::stod(argv[5]);
    if(argc > 6) h        = std::stod(argv[6]);
    if(argc > 7) experiment_id = std::stoi(argv[7]);

    // auto sites = SpinHalf(N,{"ConserveQNs=",false});

    // 初态：|00...0> 的 MPS
    // InitState init(sites);
    // for(int i = 1; i <= N; ++i) init.set(i,"Dn");
    // MPS psi(init);

    // 构造 GHZ
    // makeGHZMPS(psi,sites,maxDim);

    double J = 1.0;
    // double h = 1.0;
    MPS psi;
    std::ostringstream oss;
    oss << "MPS/ising_N" << N
        << "_J" << std::fixed << std::setprecision(2) << J
        << "_h" << std::fixed << std::setprecision(2) << h
        << ".mps";

    std::string fname = oss.str();
    readFromFile(fname, psi);
    if(length(psi) != N)
    {
        std::cerr << "Error: loaded psi length = " << length(psi)
                << " but N = " << N << "\n";
        return 1;
    }
    // MPS psi = makeTFIMGroundState(sites, J, h);
    // brickwork 浅层 Haar 随机线路
    std::random_device rd;
    std::mt19937 rng(rd());
    applyBrickworkCircuit(psi,N,depth,rng,maxDim);

    // 一次性算出整个标准基概率分布
    auto probs = computeProbDistribution(psi);
    std::size_t dim = (std::size_t(1) << N);
    applyDepolarizationProb(probs, dprob);
    // ===== 统计每个结果出现次数 =====
    std::vector<std::size_t> counts(dim,0);

    // 使用同一个离散分布对象做多次采样
    std::discrete_distribution<std::size_t> dist(probs.begin(),probs.end());

std::vector<std::size_t> checklist = {10000, 20000, 50000, 100000, 200000, 500000, 1000000};
std::size_t check_idx = 0; // 指向下一个要触发的 checklist 元素

auto write_checkpoint = [&](std::size_t curSamples) -> int {
    // === 计算 sum_C3 ===
    unsigned long long sum_C3 = 0ULL;
    unsigned long long sum_C2 = 0ULL;
    for (std::size_t x = 0; x < dim; ++x) {
        sum_C3 += choose3(counts[x]);
        sum_C2 += choose2(counts[x]);
    }

    long double S_ld = (long double)sum_C3;
    long double S_nd = (long double)sum_C2;
    // d = 2^N
    long double d = powl(2.0L, (long double)N);

    // numerator = S * (d+2)(d+1)
    long double numerator = S_ld * (d + 2.0L) * (d + 1.0L);
    long double numerator_nd = S_nd * (d + 1.0L);

    // denominator = C(curSamples,3)
    long double sample_ld = (long double)curSamples;
    long double denom = sample_ld * (sample_ld - 1.0L) * (sample_ld - 2.0L) / 6.0L;
    long double denom_nd = sample_ld * (sample_ld - 1.0L) / 2.0L;

    long double F = numerator / denom;
    long double Purity = numerator_nd / denom_nd;

    

    // === 写入文件：nSamples 用 curSamples（也就是 s+1）===
    std::ostringstream oss;
    oss << "TIM_complete/n" << N
        << "_depth" << depth
        << "_J" << std::fixed << std::setprecision(2) << J
        << "_h" << std::fixed << std::setprecision(2) << h
        << "_p" << std::fixed << std::setprecision(2) << dprob
        << "_Sample" << curSamples
        << ".jsonl";

    std::string filename = oss.str();

    std::ofstream ofs(filename, std::ios::app);
    if (!ofs) {
        std::cerr << "Error: cannot open " << filename << " for writing\n";
        return 1;
    }

    ofs << "{";
    ofs << "\"experiment_id\":" << experiment_id << ",";
    ofs << "\"N\":" << N << ",";
    ofs << "\"depth\":" << depth << ",";
    ofs << "\"nSamples\":" << curSamples << ",";
    ofs << "\"MaxDim\":" << maxDim << ",";
    ofs << "\"p\":" << std::setprecision(2) << dprob << ",";
    ofs << "\"Purity\":" << std::setprecision(20) << Purity << ",";
    ofs << "\"F\":" << std::setprecision(20) << F;
    ofs << "}\n";

    return 0;
};

    for(int s = 0; s < nSamples; ++s)
    {
        std::size_t x = dist(rng); // x 是 0..2^N-1 的整数，对应某个 bitstring
        // 统计次数
        counts[x]++;
        std::size_t curSamples = s + 1;
        if (check_idx < checklist.size() && curSamples == checklist[check_idx]) {
            if (write_checkpoint(curSamples) != 0) return 1;
            ++check_idx;
        }
        
    }


    return 0;
}
