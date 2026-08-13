# Causal audit--return memory benchmark v0.1

Status: preregistration-ready theoretical benchmark and simulator. It has not
yet been demonstrated on physical hardware.

## 1. What the benchmark measures

The benchmark asks whether a sequential device preserves a coherent quantum
state between temporally separated interactions in a way that unlimited
adaptive classical memory cannot reproduce.

For each of `n` slots, a trusted verifier prepares

\[
 |\Phi^+\rangle_{R_iA_i}
 =\frac{|00\rangle+|11\rangle}{\sqrt2},
\]

sends `A_i` to the device, receives `B_i`, and immediately sequesters `B_i`.
The next input is not released before the previous output is committed. After
all slots, the device also commits its terminal memory port to the verifier;
this may be a coherent system in the unrestricted experiment but is restricted
to an orthogonal classical transcript under the null. A hidden random choice
then selects one branch:

1. **AUDIT.** Measure every reference `R_i` in `Z`. Decode a parity guess from
   the committed terminal memory port. Sequestered outputs are unavailable.
   Under the null this is exactly a decoder of the classical transcript.
2. **RETURN.** Apply the preregistered joint decoder to the terminal memory
   port and sequestered outputs, then test
   `(|Phi+><Phi+|)^tensor n` and a standard visible reset state for the terminal
   port. The reset does not certify deletion of an inaccessible purification.

The measured probabilities are `P_A` and `F_R`. With a preregistered weight
`lambda`, the score is

\[
 S_\lambda=\lambda P_A+(1-\lambda)F_R.
\]

For the default balanced benchmark, `lambda=1/2`.

## 2. Exact classical-memory null

The null permits:

- arbitrary local quantum instruments in each slot;
- arbitrary finite classical outcomes;
- disposable within-slot quantum ancillas;
- unlimited adaptive classical memory and shared classical randomness; and
- arbitrary transcript-conditioned joint recovery in RETURN.

It forbids:

- any coherent quantum state persisting from one slot to the next;
- entanglement shared among slot ancillas before the stream;
- later access to an already sequestered output; and
- access to returned carriers in AUDIT.

Under these assumptions the exact null support is

\[
 C_{n,\lambda}
 =\max_{0\le t\le1}
 \left[
 \frac{\lambda}{2}(1+t^n)
 +(1-\lambda)
 \left(\frac{1+\sqrt{1-t^2}}2\right)^n
 \right].
\]

For every `n>=2` at `lambda=1/2`, `C=0.75`. This is the primary certification
threshold. It is not a generic Bell inequality or a device-independent bound;
it is a theorem for the declared causal interface.

For comparison, if all carriers can be processed jointly before committing a
classical outcome, the bound becomes

\[
 C_{\mathrm{collective}}(\lambda)
 =\frac12+\frac12\sqrt{\lambda^2+(1-\lambda)^2},
\]

which is approximately `0.853553` at equal weight. A coherent device with one
persistent parity-accumulator qubit can attain `(P_A,F_R)=(1,1)`.

## 3. Estimation and certification rule

Let `N_A,N_R` be fixed sample sizes and `K_A,K_R` the accepted AUDIT and RETURN
trials. Define

\[
 \widehat S
 =\lambda\frac{K_A}{N_A}
 +(1-\lambda)\frac{K_R}{N_R}.
\]

For independent Bernoulli trials, a one-sided weighted Hoeffding radius is

\[
 r_\alpha=\sqrt{\frac12
 \left(\frac{\lambda^2}{N_A}
 +\frac{(1-\lambda)^2}{N_R}\right)
 \log\frac1\alpha}.
\]

Before seeing outcomes, preregister:

- upward-bias bounds `delta_A` and `delta_R` for the two estimators;
- a null-model allowance `delta_null` for source and isolation uncertainty;
- `n`, `lambda`, `alpha`, sample sizes, decoder, exclusions, and stopping rule.

Certify only if

\[
 \widehat S-\lambda\delta_A-(1-\lambda)\delta_R-r_\alpha
 > C_{n,\lambda}+\delta_{\mathrm{null}}.
\]

Under the fixed-sample independent-trial model, this has false-positive
probability at most `alpha`. Optional stopping, choosing `lambda` after seeing
data, or choosing among several lengths requires a sequential-valid or
multiple-testing correction not supplied by this version.

## 4. Power and shot planning

For a declared alternative score `S_alt`, define the adjusted gap

\[
 g=S_{\mathrm{alt}}
 -\lambda\delta_A-(1-\lambda)\delta_R
 -C_{n,\lambda}-\delta_{\mathrm{null}}.
\]

If `g<=0`, no amount of data can establish the claim under those allowances.
For `g>0`, the bundled planner chooses fixed sample sizes satisfying

\[
 r_\alpha+r_\beta<g,
\]

which guarantees power at least `1-beta` by Hoeffding. Allocation approaches
`N_A:N_R=lambda:(1-lambda)`.

Examples:

```bash
python scripts/carmenq_cli.py bound --steps 8

python scripts/carmenq_cli.py plan \
  --steps 8 --audit-probability 0.97 --return-fidelity 0.95 \
  --alpha 0.01 --beta 0.1 \
  --audit-systematic 0.005 --return-systematic 0.005 \
  --null-slack 0.005

python scripts/carmenq_cli.py analyse \
  --steps 8 --audit-successes 9700 --audit-trials 10000 \
  --return-successes 9500 --return-trials 10000 \
  --alpha 0.01 --audit-systematic 0.005 \
  --return-systematic 0.005 --null-slack 0.005
```

## 5. Required controls

At minimum:

1. **No coherent memory:** remove or dephase the persistent memory carrier.
2. **Strong classical record:** measure every event and retain the transcript;
   AUDIT should improve while global RETURN collapses.
3. **Collective-access control:** deliberately withhold sequestration to verify
   the higher collective classical-record curve.
4. **Output-revisit control:** allow one earlier carrier to return to the
   device; this demonstrates why immediate sequestration is causal, not
   cosmetic.
5. **Decoder-null control:** replace the coherent inverse with identity or a
   mismatched decoder.
6. **Reference-source control:** characterize EPR preparation and measurement
   bias independently.
7. **Blind randomized branch choice:** the device must commit before learning
   which branch will be scored.
8. **Leakage audit:** bound unintended coherent communication between slots,
   including buses, clocks, resonators, shared ancillas, and pulse tails.

## 6. Hardware realization

The minimal honest device has one persistent memory qubit `M`. In slot `i`, it
applies `CNOT(A_i -> M)` and returns `A_i`; after the last slot it hands `M` to
the verifier. AUDIT measures `M` in `Z`. RETURN applies `CNOT(B_i -> M)` for
every sequestered output and then tests all EPR pairs plus the visible reset
`M=0`.

Candidate platforms include superconducting circuits with a bus/memory mode,
trapped ions with one long-lived accumulator ion, neutral-atom registers, and
photonic time bins coupled to a matter memory. Platform selection must account
for whether outputs can actually be isolated and revisited by the trusted
recovery station without returning them to the tested device.

## 7. Forecast model and functional length

The bundled forecast uses the explicit replaceable law

\[
 P_A(n)=\frac{1+c_0c^n}{2},\qquad
 F_R(n)=f_0f^n.
\]

Its defaults (`c0=0.995`, `c=0.998`, `f0=0.995`, `f=0.997`) are illustrative,
not calibrated hardware specifications. The generated data show expected
score, classical margins, robust adjusted margins, and conservative required
shots for `n=1,...,200`. Under the illustrative defaults, the raw expected
score first falls below the balanced streaming null at `n=151`; with the
documented `0.005` per-branch bias allowance and `0.005` null slack,
fixed-sample planning ceases to be feasible at `n=144`. These numbers describe
the example law, not any hardware platform.

Define the **certifiable functional coherent-memory length** as the largest
preregistered `n` for which the lower confidence score exceeds the enlarged
classical null. This is more operational than a bare `T2`: it requires the
memory to retain a temporal predicate and subsequently participate in coherent
global return.

## 8. Interpretation limits

A positive result certifies incompatibility with the declared classical-memory
comb under trusted source, sequestration, and measurement assumptions. It does
not identify a unique microscopic implementation, prove consciousness, favor
Everett, establish global erasure of every environment, or certify an
unrestricted hidden device. A failed result may reflect ordinary noise or
insufficient statistical power rather than absence of quantum memory.
