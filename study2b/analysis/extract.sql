-- Sense Foraging Study 2b: every fact the analysis needs, one row per fact. READ-ONLY.
-- Run through the Supabase CLI from the radlab checkout (README):
--   npx supabase db query --linked -o json -f <path>/extract.sql > F:/radlab_data/sf_study2b/raw/extract_<date>.json
--
-- `source` says what a row is:
--   enrollment     one per enrollment and schedule row (consent, session status, step crashes)
--   questionnaire  every submission (questionnaire_responses), with step_index (the presented
--                  position in the session: 0 background, 1 vignettes, 2 the scale, 3 the trial
--                  statements, 4-12 the first half of the validity block in a per-participant
--                  order, 13 the pause offer, 14-19 the second half, 20 debrief), received_at
--                  (server clock), resubmission_of (the platform's flag), presentation_format,
--                  and the answers (which carry _timing and, for the vignettes, _order; the
--                  pause offer's answer is its whole record, pause_1)
--   vas            the life-satisfaction face scale: an activity, so it lands in vas_responses
--   step           each exited step (participant_step_timings): entry and exit times
-- `participant` is the random code (study_enrollments.external_id). Prolific IDs live only
-- in external_identities, which this query never touches. Staff test enrollments are
-- included, marked is_test, so 01_prepare can count their exclusion.
WITH es AS (
  SELECT e.external_id                                                 AS participant,
         (coalesce(e.is_test, false) OR e.external_id LIKE 'test-%')   AS is_test,
         e.external_meta->>'prolific_study_id'                         AS prolific_study_id,
         e.repository_consent, e.status AS enrollment_status, e.consent_date, e.enrolled_at,
         ps.id AS schedule_uuid, ps.status AS session_status, ps.completed_at AS session_completed_at,
         (SELECT count(*) FROM session_diagnostics d WHERE d.schedule_id = ps.id AND d.kind = 'step_crash') AS step_crashes
    FROM study_enrollments e
    LEFT JOIN participant_schedule ps ON ps.participant_id = e.profile_id AND ps.study_id = e.study_id
   WHERE e.study_id = '1864d663-d01c-5d82-91f0-3047462470c3'
)
SELECT 'enrollment' AS source, es.participant, es.is_test, es.prolific_study_id, es.repository_consent,
       es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text AS schedule_id, es.session_status, es.session_completed_at,
       NULL::text AS response_id, NULL::text AS slug, NULL::int AS step_index, NULL::timestamptz AS received_at,
       NULL::text AS resubmission_of, NULL::text AS presentation_format, NULL::jsonb AS responses, NULL::int AS value,
       NULL::timestamptz AS entered_at, NULL::timestamptz AS exited_at, es.step_crashes
  FROM es
UNION ALL
SELECT 'questionnaire', es.participant, es.is_test, es.prolific_study_id, es.repository_consent,
       es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text, es.session_status, es.session_completed_at,
       r.id::text, r.questionnaire_slug, r.step_index, r.received_at,
       r.resubmission_of::text, r.presentation_format, r.responses, NULL::int,
       NULL::timestamptz, NULL::timestamptz, es.step_crashes
  FROM es
  JOIN questionnaire_responses r ON r.schedule_id = es.schedule_uuid
UNION ALL
SELECT 'vas', es.participant, es.is_test, es.prolific_study_id, es.repository_consent,
       es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text, es.session_status, es.session_completed_at,
       v.id::text, vs.slug, v.step_index, coalesce(v.received_at, v.responded_at),
       v.resubmission_of::text, NULL::text, NULL::jsonb, v.value::int,
       NULL::timestamptz, NULL::timestamptz, es.step_crashes
  FROM es
  JOIN vas_responses v ON v.schedule_id = es.schedule_uuid
  JOIN vas_scales vs ON vs.id = v.scale_id
UNION ALL
SELECT 'step', es.participant, es.is_test, es.prolific_study_id, es.repository_consent,
       es.enrollment_status, es.consent_date, es.enrolled_at,
       es.schedule_uuid::text, es.session_status, es.session_completed_at,
       NULL::text, t.subcategory, t.step_index, NULL::timestamptz,
       NULL::text, NULL::text, NULL::jsonb, NULL::int,
       t.entered_at, t.exited_at, es.step_crashes
  FROM es
  JOIN participant_step_timings t ON t.participant_schedule_id = es.schedule_uuid
 ORDER BY 2, 1, 14;
