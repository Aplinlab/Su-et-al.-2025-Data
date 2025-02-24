"""Functions specific to the analysis of amplitude sweeps."""

import base
import classes
import constants

def max_amp(stim_type, markers, p):
    max_amp = 1
    for marker in markers:
        m = p.search(marker[0])
        max_amp = max(max_amp, round(float(m.group(f"{stim_type}val")), 1))
    return max_amp


def analyse_amplsweep(filename, tick_dt, data_signal, data_mechstim, data_elecstim, markers, p, response_window_ms, threshold_height):
    df_summary = base.output_table("freq")

    correction_factor_mech = max_amp("mech", markers, p)
    correction_factor_elec = max_amp("elec", markers, p)
    pkval_mech = base.trigger_value(data_mechstim, correction_factor_mech)
    pkval_elec = base.trigger_value(data_elecstim, correction_factor_elec)

    for i,marker in enumerate(markers):
        start = int(marker[1]/tick_dt)
        if i < len(markers)-1:
            end = int(markers[i+1][1]/tick_dt)
        else:
            end = len(data_signal)

        mech_val, elec_val = base.test_values(filename, marker[0], p, is_amplsweep=True)
        if max(mech_val, elec_val) > 2:
            raise Exception(f"Amplitude has been set higher than 2 ({marker[0]} in {filename}). This interferes with trigger detection.")
        elif max(mech_val, elec_val) < 0.5:
            raise Exception(f"Amplitude has been set lower than 0.5 ({marker[0]} in {filename}). This interferes with trigger detection.")
        elif min(mech_val, elec_val) != 0:
            raise Exception(f"Neither amplitude value is 0 ({marker[0]} in {filename}). This interferes with trigger detection.")
        trigger_height_adjustment = min(1, max(mech_val, elec_val))
        triggers_mech = base.triggers(data_mechstim[start:end], pkval_mech*0.9*trigger_height_adjustment, 7)
        triggers_elec = base.triggers(data_elecstim[start:end], pkval_elec*0.9*trigger_height_adjustment, 3)
        stim_type, stim_value, triggers_phases = base.separate_sweep_phases("ampl", filename, mech_val, elec_val, triggers_mech, triggers_elec)

        responses_condition = base.detect_responses(data_signal[start:end], triggers_phases[0], response_window_ms, tick_dt, threshold_height)
        responses_mech_itlv = base.detect_responses(data_signal[start:end], triggers_phases[1], response_window_ms, tick_dt, threshold_height)
        responses_elec_itlv = base.detect_responses(data_signal[start:end], triggers_phases[2], response_window_ms, tick_dt, threshold_height)
        responses_mech_rcvr = base.detect_responses(data_signal[start:end], triggers_phases[3], response_window_ms, tick_dt, threshold_height)
        responses_elec_rcvr = base.detect_responses(data_signal[start:end], triggers_phases[4], response_window_ms, tick_dt, threshold_height)

        # TODO
        # for each trigger, look at data within a set window
        # run peak detect using a threshold multiple of SD
        # from that, calculate rate, first response, and first 10 responses

        # return "amplitude", stim_type, stim_value, rate_cond, rate_itlv_mech, rate_itlv_elec, first_itlv_mech, first_itlv_elec, firstten_itlv_mech, firstten_itlv_elec, firstten_rcvr_mech, firstten_rcvr_elec

    return df_summary
