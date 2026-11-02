# Study 2b: the free-choice pause (behavioural criterion)

Decision 10 in `decisions_after_2a.md`. Built 2 November 2026 as the platform component `sensing_pause` (radlab repo, `src/components/questionnaire/composable/SensingPauseOffer.jsx`, commit 64a87ccf on `dev`). Preview, nothing saved: radlab.zone dev site, `/dev/sensing-pause-preview` (`?fast=1` shortens the timers).

## How the design was chosen

Five designs were written up and judged by an independent reviewer (a separate agent given the designs and the criteria, not a preference):

| Rank | Design | Verdict |
|---|---|---|
| 1 | Repeated offers | Only design with a reliability estimate and a way to separate fatigue breaks from sensing; cut from three offers to two |
| 2 | One offer, timed pause | Right construct (choosing non-goal time at a real cost on fixed pay), but one binary trial is too unreliable for a confirmatory test |
| 3 | Waiting room (forced wait, choose sensing or reading) | Holds time cost constant, but a forced wait annoys, borders on deception, and depends on how interesting the article is |
| 4 | Sense menu (pick sensory invitations) | Prompts turn sensing into a task; measures compliance and completionism |
| 5 | Hold to stay | Holding a button is a motor task that defeats the construct; poor accessibility |

## The design

Two offers, each one survey page:
- **Offer 1** at about a third of the way through, before the scale, so its statements cannot prime acceptance. **[Norm: before or after the vignettes? The vignettes name the two modes.]**
- **Offer 2** at about three-quarters, mid validity block, at least 5 minutes after the scale.

Each offer: "If you'd like, you can take a short break here to look, listen and feel what's around you. There's nothing to do and nothing to report. It doesn't affect your payment." Two equal bordered buttons, **Take a pause** and **Continue the survey**; which side each is on is random per participant and the same for both offers.

- **Take a pause:** a quiet screen with one line ("Take as long as you like. Tap 'I'm back' when you're ready."), a slow colour field that cues nothing (no breathing rhythm, no focus point; still under reduced-motion settings) and **I'm back**. After 3 minutes on the page the line becomes "Whenever you're ready." After 10 minutes in all the pause ends itself (`ended_by = "cap"`). Then "Anything you noticed? (optional)".
- **Continue the survey:** "Any reason? (optional)": No time right now / Not a good place / Not interested / Something else / Skip.
- **After offer 2 only:** "Where are you right now?" (home / work or school / outdoors / travelling / somewhere else) and "How tired are you right now?" (1 Not at all to 6 Very), each with Prefer not to answer. Asked after the second decision, so they cannot steer either decision.

## What is recorded (per offer)

`choice` (pause / continue), `decision_ms`, `decline_reason`, `note`, `pause_button_side`, `device` (phone / tablet / desktop), `viewport_w`, `local_hour`; for a pause, `pause.total_ms` split into `on_page_ms` (tab visible and focused), `unfocused_ms` (visible, another window in front) and `hidden_ms` (another tab or app, phone locked), with episode counts, `taps`, `keys`, `scrolls`, `ended_by`, `nudged`; on offer 2, `setting` and `tired`. Time is on the browser's monotonic clock. Tested 2 November in headless Chrome: a 4-second tab switch during a pause recorded 4.0 s hidden and 4.4 s on the page.

## Criterion variables (proposed for the preregistration)

- **Primary: engaged-pause count, 0–2.** A pause counts when it was taken, lasted at least 20 s on the page, and time away (unfocused + hidden) was under half the total. Tested against the scale score by Spearman's ρ, with an ordinal regression controlling for device.
- **Secondary:** acceptance at offer 1 alone (0/1), the decision least affected by fatigue; total on-page seconds, log(1 + x), with "not interested" decliners compared against "no time / not a good place" decliners (the scale should predict "won't" more than "can't").
- **Discriminant check:** the criterion relates more to the scale than to tiredness.
- **Exploratory:** sensory detail in the notes, coded blind to scale score.

## Risks

1. **Too few (or too many) pause.** Two offers raise the chance of any pause. Proposed: a pilot of about 60, and if offer 1 acceptance is below 10% or above 90%, the wording is revised once before the main launch; the fallback analysis is on-page seconds.
2. **Pausing for reasons other than sensing** (tiredness, a break, the phone). The on-page / away split removes walking away; offer 1 comes before fatigue and the scale; tiredness, setting and decline reasons are recorded.

## Still to do before launch

- Export: flatten `sensing_pause` answers into columns (`pause_1_choice`, `pause_1_on_page_s`, …) in the study export.
- Ethics amendment: describe the offers and the recorded timing.
- Placement and the 20 s / 50% thresholds into the preregistration; the pilot rule, if adopted.
