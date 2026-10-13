-- Sense Foraging Study 2 and its two-week retest: every fact the analysis needs, one
-- row per fact. READ-ONLY. Run through the Supabase CLI (README):
--   npx supabase db query --linked -o json -f <path>/extract.sql > F:/radlab_data/sf_study2/raw/extract_<date>.json
--
-- `source` says what a row is:
--   enrollment     one per enrollment and schedule row (posting, consent, session status)
--   questionnaire  every submission (questionnaire_responses), with step_index (the presented
--                  position in the session: 0 background, 1 short form, 2 candidates, 3-13 the
--                  shuffled validity block, 14 debrief), received_at (server clock),
--                  resubmission_of (the platform's flag), presentation_format, and the answers
--                  (which carry _timing and, for the candidates, _order.positions)
--   vas            the life-satisfaction face scale: an activity, so it lands in vas_responses
--   step           each exited step (participant_step_timings): entry and exit times
-- `study` is 'main' (Study 2) or 'retest'. `participant` is the Study 2 random code
-- (study_enrollments.external_id). For retest rows it is the LINKED Study 2 code, from
-- external_identity_links, which deidentify_external_enrollments writes (in either
-- direction) just before it deletes the Prolific IDs; until it has run, retest rows are
-- unlinked (participant NULL, n_links 0). Prolific IDs (external_identities) are never read.
-- Staff test enrollments are included, marked is_test, so 01_prepare can count their exclusion.
WITH ids AS (
  SELECT '528674d9-f691-5cff-a81c-1f4f0996b08b'::uuid AS main_id,     -- Sense Foraging Study 2
         '65ffadea-be05-56c7-922a-20b78258c6ab'::uuid AS retest_id    -- Sense Foraging Study 2 retest
),
links AS (
  SELECT l.surrogate_a AS main_code, l.surrogate_b AS retest_code
    FROM external_identity_links l, ids
   WHERE l.study_a = ids.main_id AND l.study_b = ids.retest_id
  UNION
  SELECT l.surrogate_b AS main_code, l.surrogate_a AS retest_code
    FROM external_identity_links l, ids
   WHERE l.study_a = ids.retest_id AND l.study_b = ids.main_id
),
enr AS (
  SELECT e.id AS enrollment_id, e.profile_id, e.study_id,
         CASE WHEN e.study_id = ids.main_id THEN 'main' ELSE 'retest' END            AS study,
         CASE WHEN e.study_id = ids.main_id THEN e.external_id
              ELSE (SELECT min(k.main_code) FROM links k WHERE k.retest_code = e.external_id) END AS participant,
         CASE WHEN e.study_id = ids.retest_id THEN e.external_id END                  AS retest_code,
         CASE WHEN e.study_id = ids.retest_id
              THEN (SELECT count(*) FROM links k WHERE k.retest_code = e.external_id) END AS n_links,
         (coalesce(e.is_test, false) OR e.external_id LIKE 'test-%')                  AS is_test,
         e.external_meta->>'prolific_study_id'                                        AS prolific_study_id,
         e.repository_consent,
         e.status                                                                     AS enrollment_status,
         e.consent_date,
         e.enrolled_at
    FROM study_enrollments e, ids
   WHERE e.study_id IN (ids.main_id, ids.retest_id)
),
es AS (
  SELECT enr.*, ps.id AS schedule_uuid, ps.status AS session_status, ps.completed_at AS session_completed_at,
         (SELECT count(*) FROM session_diagnostics d WHERE d.schedule_id = ps.id AND d.kind = 'step_crash') AS step_crashes
    FROM enr
    LEFT JOIN participant_schedule ps ON ps.participant_id = enr.profile_id AND ps.study_id = enr.study_id
)
SELECT 'enrollment' AS source, es.study, es.participant, es.retest_code, es.n_links, es.is_test, es.prolific_study_id,
       es.repository_consent, es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text AS schedule_id, es.session_status, es.session_completed_at,
       NULL::text AS response_id, NULL::text AS slug, NULL::int AS step_index, NULL::timestamptz AS received_at,
       NULL::text AS resubmission_of, NULL::text AS presentation_format, NULL::jsonb AS responses, NULL::int AS value,
       NULL::timestamptz AS entered_at, NULL::timestamptz AS exited_at, es.step_crashes
  FROM es
UNION ALL
SELECT 'questionnaire', es.study, es.participant, es.retest_code, es.n_links, es.is_test, es.prolific_study_id,
       es.repository_consent, es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text, es.session_status, es.session_completed_at,
       r.id::text, r.questionnaire_slug, r.step_index, r.received_at,
       r.resubmission_of::text, r.presentation_format, r.responses, NULL::int,
       NULL::timestamptz, NULL::timestamptz, es.step_crashes
  FROM es
  JOIN questionnaire_responses r ON r.schedule_id = es.schedule_uuid
UNION ALL
SELECT 'vas', es.study, es.participant, es.retest_code, es.n_links, es.is_test, es.prolific_study_id,
       es.repository_consent, es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text, es.session_status, es.session_completed_at,
       v.id::text, vs.slug, v.step_index, coalesce(v.received_at, v.responded_at),
       v.resubmission_of::text, NULL::text, NULL::jsonb, v.value::int,
       NULL::timestamptz, NULL::timestamptz, es.step_crashes
  FROM es
  JOIN vas_responses v ON v.schedule_id = es.schedule_uuid
  JOIN vas_scales vs ON vs.id = v.scale_id
UNION ALL
SELECT 'step', es.study, es.participant, es.retest_code, es.n_links, es.is_test, es.prolific_study_id,
       es.repository_consent, es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text, es.session_status, es.session_completed_at,
       NULL::text, t.subcategory, t.step_index, NULL::timestamptz,
       NULL::text, NULL::text, NULL::jsonb, NULL::int,
       t.entered_at, t.exited_at, es.step_crashes
  FROM es
  JOIN participant_step_timings t ON t.participant_schedule_id = es.schedule_uuid
 ORDER BY 2, 3, 4, 1, 18;
