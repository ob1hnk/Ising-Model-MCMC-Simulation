https://courses.physics.illinois.edu/phys498cmp/sp2022/Ising/IsingModel.html

# I. Overview

우리의 목적은 renormalization group을 이해하는 것이고, 그것은 다음의 질문을 답한다: **How do the rules of the universe change when you look at in a coarse-grain way?**

# II. Simulating an Ising Model

## 1. Minimal Background

### The Ising Model

결정(crystal) 물질의 자기적 성질을 격자 상의 자기 모먼트 방향으로 나타내는 간단한 모델이다.

$$
E = - \sum_{i,j} J_{ij} \, s_i s_j + \sum_i h_i \, s_i
$$

$s_i$: 각 모멘트의 가능한 방향.

$J_{ij}$: interaction matrix의 요소. 모멘트 i, j 사이의 상호작용 에너지.

$h_i$: 각 위치 i에서의 자기장.

이징 모델은 $s_i = \pm 1$임. 즉, 이웃과만 상호작용함.

이징 모델을 구성(configuration)한다는 건 모든 모멘트의 방향을 명시하는 것임.
이때 ‘구성’은 ‘snapshot’과도 같다. 예를 들어 $i$번째 항목에 $\pm 1$을 포함하는 벡터는 모멘트 $i$가 위를 가리키는지 아래를 가르키는지 지정함.

### Statistical Mechanics

1. 아이징 모델에서, 계의 configuration $c$ 는 각 스핀의 방향이다.
2. 각 configuration $c$ 는 에너지 $E(c)$ 가 있다.
3. 온도 $T$ 의 energy reservoir과 맞닿아있는 고전적인 계에서, configuration이 $c$ 일 확률 $P$ 는
    
    $
    P(c) = \frac{e^{-\beta E(c)}}{\sum_{c} e^{-\beta E(c)}}\ , \ \beta \equiv 1/(kT)
    $
    
    이다. $k$는 볼츠만 상수로, 우리는 $k = 1$만 다룰 것이다.
    

### Markov Chains

위에서 봤듯이, 
$𝑃(𝑐)∼exp(−𝛽𝐸(𝑐))$ 이다. 먼저, $c$ 를 찾는 알고리즘 중 하나인 Markov Chain Monte Carlo 알고리즘을 짜자.

Markov chain은 상태 $c$ (예: configuration)에 있을 때 $c$에 의존하는 특정 확률로 다른 상태 $c'$로 가는 것을 나타낸다. $P(c →𝑐′)$은 $c$ 에 의존함. 이를 memoryless process라고 부른다.

memoryless한 방식으로 (non-pathalogical) random walk를 충분히 많이 하면 반드시 
$𝜋(𝑐)$의 확률로 상태 $c$에 도달한다. (
$𝜋$: stationary distribution) 서로 다른 markov chain은 서로 다른 
$𝜋$ 값을 가진다.

Markov chain이 non-pathological 하려면:

1. *aperiodic*; it doesn’t cycle between configurations in a subset of the system’s full configuration space.
2. *connected*; given two configurations, there is a path with non-zero probability that the chain *could* (but not necessarily *will*) follow to get from one to the other.

아이징 모델을 simulate하기 위해,  $\pi(c) \sim e^{-\beta Ec}$ 와 같은 stationary distribution을 가지는 markov chain을 만들자. 그때 Metropolis-Hastings algorithm을 이용한다.

### The Metropolis-Hastings Algorithm

우리의 desired stationary distribution을 알고 있을 때, Metropolis-Hastings algorithm을 써서 markov chain을 만든다.

**방법**

1. 어떤 configuration $c$ 에서 시작
2. $T(c →𝑐′)$의 확률로 새로운 trial configuration $c'$으로 가는 이동을 설정
3. $\min\left( 1,\; \frac{\pi(c')\, T(c' \to c)}{\pi(c)\, T(c \to c')} \right)$의 확률로 $c'$으로 이동하기

2.1.  $T(c →𝑐′)$를 선택하는 방법

이동의 과정을 선택하는 것은 MCMC의 핵심이라고 할 수 있다. 각자 장단점이 있는 선택지가 여러개 있음. 우리의 경우 가는 확률과 오는 확률이 동일한 것으로 설정하여, 알고리즘 구현을 단순하게 할 수 있도록 함. $𝑇(𝑐→𝑐′)=𝑇(𝑐′→𝑐)$

3.1. Acceptance condition

원하는 분포 𝜋에 점근하는 markov chain을 configure하려면 알고리즘은 configuration space을 탐색하는 과정에서 𝜋을 어떤 방식으로든 포함해야 한다. 메트로폴리스 알고리즘은 acceptance condition에서 포함함.

$$
\alpha≡\frac{\pi(c')\, T(c' \to c)}{\pi(c)\, T(c \to c')}
$$

𝛼 > 1: $c'$ is more likely

𝛼 < 1: $c$ is more likely

## 2. States and Parameters

Ising Model의 에너지 공식 revisited (below):

$$
E = - \sum_{i,j} J_{ij} \, s_i s_j + \sum_i h_i \, s_i
$$

여기에서, 우리가 다룰 케이스는 아래와 같다.

$h_i = 0$, 자기장 없음

$J_{ij} = J$ for 이웃하는 i, j. 이웃하지 않는 경우 $J_{ij} = 0$

따라서 에너지 공식은 다음과 같다.

$$
E = -J \sum_{(x,y)} s_{(x,y)}\, s_{(x+a, y+b)}
, \ where \ (𝑎,𝑏)∈{(1,0),(0,1),(−1,0),(0,−1)}
$$

우리는 $L \times L$의 그리드( $L = 3, 27, 81$ 등)를 사용할 것이다.

참고: $L \times L$ 그리드의 스핀(-1 or 1)을 랜덤 초기화하기 위해 `np.random.choice((L,L), [1,-1])`을 사용한다.

## 3. Computing the Energy

1. configuration $c$의 에너지를 계산하는 함수 (예를 들면 `def Energy(spins)` )를 작성할 것.
2. 하나의 spin flip 차이를 가지는 두개의 configurations의 에너지 차이를 계산하는 함수 `def deltaE(spin_to_flip)` 를 작성할 것.
    1. 필요한 이유: 우리가 사용할 acceptance criterion은 다음의 비율이기 때문이다.
        
        $\frac{\exp\!\left(-\beta E_{\text{new}}\right)}{\exp\!\left(-\beta E_{\text{old}}\right)}= \exp\!\left(-\beta \Delta E\right)$
        
    2. 이때, $O(1)$으로 연산 가능하도록 한다. 특히, `Energy(spins)` 를 두 번 호출하여 $O(N)$이 되지 않도록 한다. 즉, `Energy(spins)` 호출 → spin flip → `Energy(spins)` 다시 호출하지 말도록.

## 4. Writing MCMC

MCMC에서, **sweep**(하나의 sweep는 $N$번의 **step**)을 여러번 반복한다. 약 20번의 sweep을 수행한 이후(equilibration time 지난 이후)부터 데이터를 수집한다. Markov chain에서 초기 조건을 잊는 데에 시간이 걸리기 때문이다.

### **Step의 과정**

1. 현재의 configuration $c$를 고려하여, $p(c →𝑐′)$의 확률을 가지는 새로운 configuration $c'$를 선택한다. 가장 간단한 방법은 random spin을 선택하여 flip하는 것이다.
우리의 경우 $T(c \to c') = T(c' \to c) = \frac{1}{N}, \ N\ is\ the\ number\ of\ spins$ 이고, acceptance probability $\alpha$에서 $T$가 상쇄된다.

1. 에너지 차이 $\Delta(c, c') \equiv E(c') - E(c)$를 계산한다.
2. $T(c' \to c)\, T(c \to c') \, \exp\!\left[-\beta \Delta(c, c')\right] > {ran}(0,1)$인 경우, 새로운 configuration $c'$로 전환하고 그렇지 않은 경우 기존의 $c$에 머문다.

### Sweep의 과정

하나의 sweep은 $N$번의 step을 반복하는 것이다. Sweep 한번마다 측정(snapshot)을 수행하자.

```python
for sweep in range(0,10000):
	for step in range(0,N):
		# flip spins
	# measure
```

3×3 그리드에서, 매 sweep마다 생성한 snapshot으로 histogram을 만들자. 이때 (↑,↓,↑...↑)와 같이 스핀 9개의 나열으로 나타낼 수 있으며, 이를 (101…1)와 같이 이진수로 표현하여 $0$ 과 $2^9 - 1$ 사이 정수로 나타낼 수 있다. 따라서 결론적으로 매 sweep마다 하나의 정수값을 얻을 것이다.

### 결과 그래프 도출

이 histogram을 이론적 계산값과 비교하자. Theory graph를 만들기 위해, 파이썬으로 $2^9$개의 configuration의 에너지를 모두 계산하고 $\frac{\exp\!\left(-\beta E(c)\right)}{Z}, \ Z\ is\ the\ normalization\$을 이용할 수 있다.

두 개의 그래프를 포개어 나타내자.

# III. Measuring the Ising Model

여기에서는 27×27 그리드에서 측정하는 것을 살필 것이다. 에너지 $E$와 자화 제곱 $M^2$을 구할 때, 항상 하나의 스핀에 대하여 구한다.

## 1. Minimal Background

온도 $T$에서 계의 관측량 $O(c)$을 측정할 수 있다. $O(c)$는 configuration을 받아 숫자를 출력한다. 관측량의 평균값은 다음과 같다.

$$
\sum_c P(c)O(c)
$$

우리는 다음의 관측량들을 살필 것이다.

**에너지 공식**

$$
E = -J \sum_{(x,y)} s_{(x,y)}\, s_{(x+a, y+b)}
, \ where \ (𝑎,𝑏)∈{(1,0),(0,1),(−1,0),(0,−1)}
$$

**에너지 평균값 공식**

$$
\langle E \rangle
= \sum_c P(c)\, E(c)
= \sum_c \frac{e^{-\beta E(c)}}{Z}\, E(c)
= \sum_E P(E)\, E
$$

이때, $Z = \sum_c e^{-\beta E(c)}$
이다.

**에너지 분포**

$$
P(E) = \sum_{c\ with\ energy\ E} P(c)
$$

**스핀 당 자화 제곱**

계가 얼마나 정렬되었는지 보여주는 지표. 모두 정렬되면 1, 완전히 무작위면 0에 가까움.

$$
M^2(c) = \left( \frac{1}{N} \sum_i s_i \right)^2
$$

**스핀 당 자화 제곱 평균**

$$
\langle M^2 \rangle= \sum_c P(c)\, M^2(c)= \sum_{M^2} P(M^2)\, M^2
$$

**자화 제곱 확률분포**

$$
P(M^2) = \sum_{\substack{c \\ \text{with squared magnetization } M^2}} P(c)
$$

### Phases

물질의 phase(상)은 관측량과 같은 order parameter(여기에서는 $M^2$) 으로 결정되며, 이때 하나의 phase에서 값이 $0$이고 다른 phase에서 값이 nonzero이다.

2D 아이징 모델에서는 2개의 phase가 있다.

- ferromagnet(강자성): 모든 스핀이 같은 방향, 열역학적 극한에서 $M^2 \neq 0$.
- paramagnet(상자성): 스핀이 무작위로 배열, 열역학적 극한에서 $M^2 = 0$.

### A note about errors

Monte Carlo는 상관된(correlated) 샘플을 생성한다.

평균을 계산할 때는 단순히 합을 내고 샘플 수로 나누면 되지만, 오차 분석에는 상관성을 반드시 고려해야 한다.

**Standard error $\sigma_{E}$**

$$
\sigma_E = \frac{\sigma}{\sqrt{N}}
$$

- $\sigma$: standard deviation (표준편차)
- $N$: 데이터 개수

이때 비상관 데이터의 경우 $N$을 그대로 쓰면 되지만, 상관 데이터의 경우 integrated autocorrelation time(적분 자기상관 시간) $\tau$을 반영하여 독립 샘플 수를 구해야 한다.

$$
N_{effective} = \frac{N}{\tau}
$$

**파이썬 코드 예시**

```python
import stats
(mean, variance, error, autocorrelation) = stats.Stats(myNumpyArrayOfData)
```

여기에서 `error` 가 **$\sigma_{E}$** (standard error)이다.

## 2. Measuring

### i. 데이터 수집

다양한 $\beta$ 값에 대해 Monte Carlo simulation을 실행하고 데이터를 수집한다. 

- $\beta J \in [0.0, 0.1, 0.2, ..., 1.0]$,
- $\beta = \infty$ (코드에서는 매우 큰 값을 사용하면 됨)

0.0부터 1.0까지 0.1 간격 그리고 매우 큰 값 하나에 대한 데이터에서, 무질서 → 도메인 형성 → 완전 정렬 과정을 볼 수 있다.

*이때 $\beta$는 inverse temperature이다.

$$
\beta = \frac{1}{k_BT},\ k_B: Boltzman\ constant
$$

각 온도에서 최종 스핀 상태를 `matshow()`으로 그리고, 에너지 $E$와 자화율 제곱 $M^2$ 데이터를 수집하여 평균값과 오차를 계산한다.

위의 $\beta$ 값들에 대해 다음을 구한다. 

- prototypical snapshot (`pylab.matshow()` 를 사용하여 시각화)
- $P(E)$ histogram
- value of $\langle E \rangle$ (and error bar)
- $P(M^2)$ histogram
- value of $\langle M^2 \rangle$ (and error bar)

### ii. Testing: 극한 상황에서의 코드 검증

시뮬레이션이 정확한지 확인하기 위해 이론적으로 정답을 아는 두 지점과 비교한다.

- $\beta = 0$: 에너지가 너무 높아서 스핀들이 제멋대로 요동치는 상태(상자성 paramagnetic)
    - 무한대 온도($T = \infty$)
    - 스핀이 1/2 확률로 무작위라 $M \approx 0$.
- $\beta = \infty$: 에너지가 가장 낮은 경계 상태로, 모든 스핀이 한 방향을 향함(강자성 ferromagnetic)
    - 절대 0도($T = 0$)
    - 스핀이 일렬로 정렬되어 $M \approx 1$ 또는 $-1$.

두 지점에 대해, histogram 그리고 $\langle E \rangle$과 $\langle M^2 \rangle$의 평균값을 구하여 데이터 값과 비교한다. 이론적 계산값을 구하는 것은 직접(해석적으로) 하거나 짧은 파이썬 코드를 이용할 수 있다.

### iii. 상전이 지점 찾기

물질이 자성을 잃거나 얻는 임계 온도를 찾는다.

x축을 $\beta$, y축을 $\langle M^2 \rangle$로 하는 그래프를 그린다. 이때 error bar를 포함한다. 그래프에서 $M^2$ 값이 급격하게 변하는 구간을 찾는다.

### iv. 비열(specific heat) 계산 및 분석

상전이 지점에서 에너지가 급격히 변하며 비열이 peak를 형성한다.

**비열 공식**

$$
C_v = \frac{\partial E}{\partial T}
$$

*데이터를 직접 미분하면 노이즈가 심하므로 `scipy.interpolate.UnivariateSpline` 을 이용하여 에너지 곡선을 매끄럽게 만든 후 미분하자.

x축을 $T$, y축을 $C_v$로 하는 그래프를 그리고 peak가 나타나는 온도를 찾는다.