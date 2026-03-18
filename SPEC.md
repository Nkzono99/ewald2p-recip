# SPEC.md

## 1. 概要

本仕様は、**点ソース近似による Coulomb 相互作用 \(1/r\)** に対して、  
**2軸周期・1軸開放（2P: two-periodic slab）** の **Ewald reciprocal-space 寄与**を評価する  
**独立した Fortran モジュール / ライブラリ**の仕様を定める。

本モジュールは BEM 本体や FMM 本体から独立した内部ライブラリとして実装し、  
境界要素メッシュそのものではなく、以下の配列 API のみを受け取る。

- source 位置 `src_pos(3,nsrc)`
- source 電荷 `src_q(nsrc)`
- target 位置 `target_pos(3,ntgt)`

本モジュールは **reciprocal-space 寄与のみ**を返す。  
real-space 側は別実装（例: Cartesian / Taylor FMM）で評価し、  
外側の adapter が両者を合成する。

---

## 2. スコープ

### 2.1 対象

初版で対象とする問題は以下に限定する。

- kernel は **Coulomb \(1/r\)** のみ
- source は **点ソース**のみ
- 周期境界は **2軸のみ**
- 残り1軸は **開放境界**
- reciprocal-space 寄与を **有限 Fourier mode 和**で直接評価する
- **中性セルのみ**対応する
- **標準 2P Ewald の \(k=0\) 項を含む**
- self interaction は **除外**する

### 2.2 非対象

初版では以下は対象外とする。

- 非中性セル
- 一様背景電荷
- charged wall による正則化
- 面要素の form factor
- 高次要素・形状関数
- FFT / PME / NFFT / NUFFT による高速化
- real-space 項の評価
- エネルギー評価
- 力・電場以外の物理量
- softening 付き kernel
- 3軸周期 / 1軸周期 / 非周期の統合実装

---

## 3. 設計方針

### 3.1 独立ライブラリ化

本モジュールは BEM 実装に直接依存しない。  
BEM 側は要素重心座標と要素電荷を点ソース近似として本モジュールへ渡す。

本モジュールは以下を `use` してはならない。

- `mesh_type`
- `sim_config`
- BEM 固有の要素型
- FMM 固有の plan / state 型

### 3.2 API 方針

公開 API は配列ベースとし、主に以下を提供する。

```fortran
call build_recip_plan(plan, options)
call update_recip_state(plan, state, src_pos, src_q)
call eval_recip_points(plan, state, target_pos, e, phi)
call eval_recip_point(plan, state, r, e, phi)
````

### 3.3 数値実装方針

* reciprocal-space は **direct mode sum**
* source 依存量は `update_recip_state` でモード係数へ圧縮
* target ごとの評価は `eval_recip_point(s)` で行う
* zero-mode は別ルーチンで扱う
* source = target の自己相互作用は明示的に除外する

---

## 4. 問題設定

### 4.1 幾何

周期方向を (x,y)、開放方向を (z) とする。
周期長を (L_x, L_y) とする。

source 点:
[
\mathbf{x}_j = (x_j, y_j, z_j)
]

target 点:
[
\mathbf{x} = (x, y, z)
]

面内ベクトル:
[
\boldsymbol{\rho} = (x,y), \qquad \boldsymbol{\rho}_j = (x_j,y_j)
]

### 4.2 中性条件

初版では source 全電荷
[
Q_{\mathrm{tot}} = \sum_{j=1}^{N} q_j
]
が 0 であることを必須とする。

数値的には

[
|Q_{\mathrm{tot}}| \le \texttt{neutrality_tol} \cdot \max\left(1, \sum_j |q_j|\right)
]

を満たすことを要求する。

この条件を満たさない場合、本モジュールはエラーを返す。

---

## 5. 物理仕様

## 5.1 Ewald 分解

全ポテンシャルは概念的に次の和として分解される。

[
\phi = \phi^{(\mathrm{real})} + \phi^{(\mathrm{recip})} + \phi^{(k=0)} - \phi^{(\mathrm{self})}
]

本モジュールが担当するのは

* reciprocal-space の (k \neq 0) 項
* zero-mode の (k=0) 項
* self 補正（potential を返す場合のみ）

である。

real-space 項は本モジュールの責務外とする。

## 5.2 逆格子

逆格子ベクトルを

[
\mathbf{k}_{hl} =
\left(
\frac{2\pi h}{L_x},
\frac{2\pi l}{L_y}
\right)
]

とする。
((h,l) = (0,0)) は nonzero mode から除外し、zero-mode として別に扱う。

## 5.3 Ewald パラメータ

分割パラメータ (\xi > 0) を用いる。
本モジュールは `options%xi` を必須入力とする。

---

## 6. zero-mode / self / 中性条件の仕様

### 6.1 zero-mode の仕様

初版では zero-mode policy は **1通りのみ**とする。

```text
zero_mode_policy = 'standard_2p_neutral'
```

この policy では、**標準 2P Ewald の (k=0) 項を必ず含める**。
ユーザーが `k=0` 項を無効化する API は提供しない。

理由:

* 2P slab では (k=0) 項が物理境界条件の一部である
* 実装者ごとの恣意的な省略を防ぐ
* real / reciprocal / zero-mode の責務を明確に保つ

### 6.2 中性条件の仕様

初版では neutrality mode は **1通りのみ**とする。

```text
neutrality_mode = 'require_neutral'
```

非中性セルが入力された場合は即時エラーとする。
内部で一様背景や壁電荷による正則化を自動で導入してはならない。

### 6.3 self 補正の仕様

初版では self policy は **1通りのみ**とする。

```text
self_policy = 'exclude'
```

意味は以下の通り。

#### field 評価

* source = target の自己寄与は含めない
* 追加のベクトル self 項は持たない

#### potential 評価

* source = target の自己寄与を含めない
* さらに標準 self 補正
  [
  \phi_{\mathrm{self},j} = -\frac{2\xi}{\sqrt{\pi}} q_j
  ]
  を加える

注:

* この self 補正は **source 点そのものを評価する場合**に意味を持つ
* 任意 target 点で potential を返すだけなら、通常は自己項除外で十分である
* 実装では `target_is_source` を判定できる API を別途持たせてもよい

---

## 7. 公開 API

## 7.1 モジュール名

初版の公開モジュール名は以下を推奨する。

```fortran
module ewald2p_reciprocal
```

## 7.2 公開手続き

```fortran
public :: build_recip_plan
public :: destroy_recip_plan
public :: update_recip_state
public :: destroy_recip_state
public :: eval_recip_point
public :: eval_recip_points
public :: check_neutrality
```

---

## 8. 型定義

## 8.1 `ewald2p_recip_options_type`

```fortran
type :: ewald2p_recip_options_type
    real(dp) :: lx
    real(dp) :: ly
    real(dp) :: xi

    integer(ip) :: h_max
    integer(ip) :: l_max

    real(dp) :: neutrality_tol = 1.0d-12

    logical :: return_potential = .false.
    logical :: enable_zero_mode = .true.   ! 初版では常に .true. を要求
end type
```

### フィールド仕様

* `lx`, `ly`

  * 周期長
  * 正でなければならない

* `xi`

  * Ewald 分割パラメータ
  * 正でなければならない

* `h_max`, `l_max`

  * direct mode sum の打ち切り
  * (h \in [-h_{\max}, h_{\max}]), (l \in [-l_{\max}, l_{\max}])
  * ((0,0)) は除外

* `neutrality_tol`

  * 中性判定用の許容値

* `return_potential`

  * `phi` を返すかどうか

* `enable_zero_mode`

  * 初版では `.true.` 固定
  * `.false.` が与えられた場合はエラーにしてよい

## 8.2 `ewald2p_recip_plan_type`

```fortran
type :: ewald2p_recip_plan_type
    logical :: is_built = .false.

    real(dp) :: lx, ly, xi
    integer(ip) :: h_max, l_max
    real(dp) :: neutrality_tol

    integer(ip) :: nk = 0
    real(dp), allocatable :: kx(:)
    real(dp), allocatable :: ky(:)
    real(dp), allocatable :: kabs(:)

    integer(ip), allocatable :: h_list(:)
    integer(ip), allocatable :: l_list(:)
end type
```

### 意味

* reciprocal 側の source 非依存前計算を保持する
* mode 一覧、波数、絶対値などを持つ
* source 座標や電荷は保持しない

## 8.3 `ewald2p_recip_state_type`

```fortran
type :: ewald2p_recip_state_type
    logical :: is_ready = .false.

    integer(ip) :: nsrc = 0

    real(dp), allocatable :: src_pos(:, :)
    real(dp), allocatable :: src_q(:)

    complex(dp), allocatable :: s_plus(:)
    complex(dp), allocatable :: s_minus(:)
end type
```

### 意味

* source 依存状態を保持する
* `src_pos`, `src_q` は必要なら保持する
* `s_plus`, `s_minus` は mode 係数
* 正確な係数定義は内部実装に委ねる

---

## 9. 関数仕様

## 9.1 `build_recip_plan`

```fortran
subroutine build_recip_plan(plan, options, ierr, message)
```

### 入力

* `options`

### 出力

* `plan`
* `ierr`
* `message`

### 処理

* オプション妥当性検査
* mode 一覧の構築
* `kx`, `ky`, `kabs` の前計算
* `plan%is_built = .true.`

### エラー条件

* `lx <= 0`
* `ly <= 0`
* `xi <= 0`
* `h_max < 0`
* `l_max < 0`
* `enable_zero_mode = .false.`

## 9.2 `update_recip_state`

```fortran
subroutine update_recip_state(plan, state, src_pos, src_q, ierr, message)
```

### 入力

* `plan`
* `src_pos(3,nsrc)`
* `src_q(nsrc)`

### 出力

* `state`
* `ierr`
* `message`

### 処理

* source 配列サイズ検査
* 中性判定
* mode ごとの source 係数を構築
* `state%is_ready = .true.`

### エラー条件

* `plan%is_built = .false.`
* 配列サイズ不整合
* 中性条件不成立

## 9.3 `eval_recip_point`

```fortran
subroutine eval_recip_point(plan, state, r, e, phi, ierr, message)
```

### 入力

* `plan`
* `state`
* `r(3)`

### 出力

* `e(3)`
* `phi`（必要な場合）
* `ierr`
* `message`

### 処理

* (x,y) を周期 box に wrap
* (k \neq 0) の reciprocal mode 和を評価
* (k=0) 項を評価
* potential が必要なら self policy に従う

## 9.4 `eval_recip_points`

```fortran
subroutine eval_recip_points(plan, state, target_pos, e, phi, ierr, message)
```

### 入力

* `target_pos(3,ntgt)`

### 出力

* `e(3,ntgt)`
* `phi(ntgt)`（必要な場合）

### 処理

* 各 target に対して `eval_recip_point` 相当の処理を行う
* OpenMP 並列を許可してよい

---

## 10. アルゴリズム

## 10.1 `build_recip_plan`

```text
1. validate options
2. enumerate (h,l) pairs in [-h_max, h_max] x [-l_max, l_max]
3. remove (0,0)
4. compute kx, ky, |k|
5. store plan
```

## 10.2 `update_recip_state`

```text
1. validate source array sizes
2. check neutrality
3. copy src_pos, src_q if needed
4. for each nonzero mode k:
     accumulate source mode coefficients
5. mark state ready
```

## 10.3 `eval_recip_point`

```text
1. wrap x,y into [0,Lx) x [0,Ly)
2. evaluate nonzero Fourier modes
3. evaluate k=0 contribution
4. if potential requested:
     apply self policy if target corresponds to source point
```

---

## 11. 数学仕様

### 11.1 nonzero mode

内部実装は、非零波数 (\mathbf{k}\neq 0) の reciprocal-space 和を評価する。
具体的な式は 2P Ewald の標準式に従う。

内部表現として、各 mode に対して source 依存の複素係数を構成し、
target 評価では

[
e^{i\mathbf{k}\cdot\boldsymbol{\rho}}
]

因子と (z) 依存関数を用いて寄与を足し合わせる。

### 11.2 zero mode

(k=0) は nonzero mode 和に含めず、**必ず別処理**とする。

### 11.3 self term

potential を source 点で評価する場合、self term は

[
\phi_{\mathrm{self}} = -\frac{2\xi}{\sqrt{\pi}}q
]

を用いる。

---

## 12. 周期 wrap の仕様

* (x) は `[0, Lx)` に wrap
* (y) は `[0, Ly)` に wrap
* (z) は wrap しない

注意:

* reciprocal 側の mode 評価では (x,y) の periodic equivalence のみを使う
* `src_pos` は必ずしも box 内に事前正規化されていなくてもよい
* 内部で phase を計算するときは、数値的一貫性のため必要に応じて wrap してよい

---

## 13. エラー仕様

## 13.1 `ierr` の推奨値

```text
0   success
1   invalid option
2   plan not built
3   state not ready
4   size mismatch
5   non-neutral system
6   zero mode disabled but unsupported
7   memory allocation failure
8   internal numerical failure
```

## 13.2 `message`

* 固定長文字列または allocatable character を許す
* 開発時のデバッグに十分な内容を返す
* 正常終了時は空文字でよい

---

## 14. 並列化

* `update_recip_state` の mode ループは OpenMP 並列化可能
* `eval_recip_points` の target ループは OpenMP 並列化可能
* mode 係数集約時の race condition を避けること
* 初版ではまず逐次版を正しく作り、その後 OpenMP 化すること

---

## 15. テスト仕様

初版では最低限以下を確認する。

### 15.1 build / update / eval の基本動作

* mode 数 0 でないこと
* source 1点では中性エラーになること
* source 2点で中性判定が通ること

### 15.2 周期対称性

* target を (x \to x + L_x) しても結果が一致
* target を (y \to y + L_y) しても結果が一致

### 15.3 reciprocal 収束

* `h_max`, `l_max` を増やしたとき解が収束すること

### 15.4 source/target 一致時

* potential を返す場合に self policy が正しく働くこと

### 15.5 real + reciprocal 合成

* 別実装の real-space と合成し、直接和または参照解と比較すること

---

## 16. ディレクトリ構成案

```text
src/
  physics/
    field_solver/
      ewald/
        ewald2p_reciprocal.f90
        ewald2p_reciprocal_build.f90
        ewald2p_reciprocal_state.f90
        ewald2p_reciprocal_eval.f90
        internal/
          ewald2p_recip_types.f90
          ewald2p_recip_kspace.f90
          ewald2p_recip_zero_mode.f90
          ewald2p_recip_wrap.f90
          ewald2p_recip_utils.f90
```

---

## 17. 実装上の禁止事項

* BEM 固有型を `use` しない
* FMM 固有型を `use` しない
* 非中性系を黙って処理しない
* zero-mode を黙って省略しない
* self 項を暗黙の約束で扱わない
* 入出力配列サイズ不整合を黙って無視しない

---

## 18. 将来拡張

以下は将来版で追加可能とする。

* FFT / PME / NFFT backend
* 非中性系の正則化
* 一様背景電荷
* charged walls
* potential energy API
* stress / virial
* panel form factor
* softening kernel 版
* mixed backend (`direct` / `nfft`)
* target/source coincidence map を用いた self 処理の厳密化

---

## 19. 初版の結論

初版モジュールは、次の1文で要約できる。

> This module evaluates the reciprocal-space contribution of the standard neutral two-periodic Ewald sum for point charges with the Coulomb kernel (1/r), including the analytic (k=0) term and excluding self-interactions. Non-neutral systems are rejected.

この方針により、real-space FMM と疎結合な独立ライブラリとして実装できる。
