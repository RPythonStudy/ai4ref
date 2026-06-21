# Data Model: KQ·PICO 정의 (Phase 1)

권위 원본 = `config/key_questions.yml` 의 `kqs:` 리스트. 각 원소 = 1 KQ 레코드.

## Entity: KQ 레코드

| 필드 | 타입 | 필수 | 설명·검증 규칙 |
|---|---|---|---|
| `kq` | str | ✅ | KQ 주제(자연어) |
| `question_type` | enum | ✅ | `intervention`\|`diagnostic`\|`prognostic`\|`predictive` |
| `design_strictness` | enum | ✗ | `strict`\|`loose`. 미지정 → 기본 `loose`(FR-007·헌법 VII) |
| `pico` | PICO | ✅ | 아래 PICO 엔티티 |
| `guideline` | Guideline | ✅ | 감시 대상 지침 앵커 |
| `guideline_refs` | list[str(PMID)] | ✗ | 검증 A: T 이전 근거(recall 참값). 없으면 채점 제외 |
| `post_guideline_landmarks` | list[str(PMID)] | ✗ | 검증 B: T 이후 랜드마크 |
| `collection` | str | ✗ | Zotero 보관 컬렉션명 |
| `enabled` | bool | ✗ | 감시 토글. 미지정 → 기본 처리(로더 정책) |
| `created_date` | str(YYYY-MM-DD) | ✗ | 생성일(메타) |

**금지 필드**: 정련된 검색식(`query` 등) 저장 금지 — 검색식은 PICO 에서 파생(헌법 III).

## Entity: PICO (OR-블록)

| 필드 | 타입 | 필수 | 설명 |
|---|---|---|---|
| `P` | list[str] | ✅ | 대상·병태 검색어 OR 리스트. 비어있으면 무효(검색식 파생 불가) |
| `I` | list[str] | ✅ | 중재/검사 검색어 OR 리스트. 비어있으면 무효 |
| `C` | str | ✗ | 대조 기준(검색식 제외 → 게이트 필터 참고) |
| `O` | str | ✗ | 결과 기준(검색식 제외 → 게이트 판정) |

- 항목 형식 = PubMed 검색어 토큰(예: `surgery[tiab]`, `"Fluid Therapy"[Mesh]`). 형식은
  파생 시 그대로 사용(이스케이프/검증은 003 검색 단계 책임).
- **정련** = `P`/`I` 에 항목 **추가**(리스트 성장 = recall 향상, FR-008).

## Entity: Guideline (지침 앵커)

| 필드 | 타입 | 필수 | 설명 |
|---|---|---|---|
| `name` | str | ✅ | 지침 이름 |
| `pmid` | str | ✗ | 지침 문헌 식별자 |
| `date` | str(YYYY-MM) | ✗ | 제정시점 T. 없으면 검증 시간축 분리(A/B) 생략 |

## 파생 규칙 (저장 안 함)

- `derived_query = build_query(pico)` = `(I[0] OR I[1] OR …) AND (P[0] OR P[1] OR …)`.
- 입력 PICO 가 고정이면 출력 결정론(재현). acceptance `recall_baseline.yml`의
  `derived_query` 와 일치해야 등가.

## 검증 규칙 요약 (로드 시, 결정론)

1. `kq`·`question_type`·`pico`·`guideline.name` 필수.
2. `pico.P`·`pico.I` 비어있지 않음.
3. 검색식 필드(`query` 등) 부재.
4. `question_type` ∈ 허용 enum, `design_strictness` ∈ {strict, loose} 또는 미지정.
5. (경고) 동일 PMID 가 `guideline_refs` ∩ `post_guideline_landmarks` 면 정의 오류 표면화.
