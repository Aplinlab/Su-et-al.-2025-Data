DISTANCE_COMPONENTS = ('sural', 'foot', 'd1', 'd2-4', 'd5')
PHASES = ('conditioning', 'interleaved', 'recovery')
SEXES = ('F', 'M')
STIMULATION_TYPES = ('electrical', 'mechanical')
TEST_CODES = {
    'long': 'long-duration',
    'freq': 'frequency',
    'ampl': 'amplitude',
    'nine': 'nine-one'
}
UNIT_TYPES = ('FA', 'SA', 'Ad', None)
UNIT_TYPES_PRINT = ('FA', 'SA', 'HT', 'IN')


#*STIMULATION PARAMETERS
PULSEWIDTH_MS = {
    'mechanical': 2.0,
    'electrical': 1.0,
    'electrical_positive': 0.2
}
DEFAULT_CONDITIONING_AMPLITUDE = 1.0
DEFAULT_CONDITIONING_FREQUENCY = 100.0
DEFAULT_CONDITIONING_FREQUENCY_MINOR = -1.0
DEFAULT_CONDITIONING_DURATION_SECONDS = 3
LONG_CONDITIONING_FREQUENCY = 50.0
LONG_CONDITIONING_DURATION_SECONDS = 300
NINEONE_CONDITIONING_FREQUENCY_MINOR = 10.0
POSTCONDITIONING_AMPLITUDE = 1.0
INTERLEAVED_DURATION_SECONDS = 3
INTERLEAVED_FREQUENCY = 25
RECOVERY_FREQUENCY = 1
SHORT_RECOVERY_DURATION_SECONDS = 20
LONG_RECOVERY_DURATION_SECONDS = 600

#*METADATA
DISTANCE_ESTIMATION_DATA = [
    {
        'age': 9,
        'sex': 'M',
        'weight': 354,
        'sural': 52,
        'foot': 22,
        'd1': 6,
        'd2-4': 19,
        'd5': 11
    },
    {
        'age': 13,
        'sex': 'M',
        'weight': 510,
        'sural': 60,
        'foot': 22,
        'd1': 8,
        'd2-4': 18,
        'd5': 11
    },
    {
        'age': 13,
        'sex': 'M',
        'weight': 558,
        'sural': 57,
        'foot': 20,
        'd1': 8,
        'd2-4': 20,
        'd5': 15
    },
    {
        'age': 14,
        'sex': 'M',
        'weight': 475,
        'sural': 60,
        'foot': 21,
        'd1': 9,
        'd2-4': 23,
        'd5': 14
    },
    {
        'age': 9,
        'sex': 'M',
        'weight': 346,
        'sural': 56,
        'foot': 20,
        'd1': 8,
        'd2-4': 21,
        'd5': 12
    },
    {
        'age': 11,
        'sex': 'M',
        'weight': 450,
        'sural': 53,
        'foot': 20,
        'd1': 8,
        'd2-4': 21,
        'd5': 13
    },
    {
        'age': 10,
        'sex': 'F',
        'weight': 259,
        'sural': 49,
        'foot': 20,
        'd1': 8,
        'd2-4': 19,
        'd5': 11
    },
    {
        'age': 14,
        'sex': 'F',
        'weight': 313,
        'sural': 54,
        'foot': 20,
        'd1': 7,
        'd2-4': 18,
        'd5': 11
    },
    {
        'age': 13,
        'sex': 'F',
        'weight': 305,
        'sural': 49,
        'foot': 22,
        'd1': 6,
        'd2-4': 16,
        'd5': 10
    },
    {
        'age': 13,
        'sex': 'F',
        'weight': 281,
        'sural': 43,
        'foot': 20,
        'd1': 6,
        'd2-4': 17,
        'd5': 10
    }
]

ANIMAL_DATA = {
    'HFF02': {
        'sex': 'M',
        'weight': 337,
        1: {
            'type': 'FA',
            'mech_amp': 2.0,
            'elec_amp': 0.8,
            'distance': {
                'sural': 1,
                'foot': 0.5
            }
        },
        2: {
            'type': 'SA',
            'mech_amp': 1.5,
            'elec_amp': 2.2,
            'distance': {
                'sural': 1
            }
        }
    },
    'HFF03': {
        'sex': 'F',
        'weight': 300,
        3: {
            'type': 'FA',
            'mech_amp': 1.5,
            'elec_amp': 1.0,
            'distance': {
                'sural': 1
            }
        }
    },
    'HFF04': {
        'sex': 'M',
        'weight': 360,
        2: {
            'type': None,
            'mech_amp': 3.0,
            'elec_amp': 2.4,
            'distance': {
                'sural': 1,
                'foot': 1,
                'd5': 1
            }
        }
    },
    'HFF05': {
        'sex': 'M',
        'weight': 369,
        1: {
            'type': None,
            'mech_amp': 1.5,
            'elec_amp': 0.6,
            'distance': {
                'sural': 1,
                'foot': 1,
                'd5': 1
            }
        },
        2: {
            'type': 'Ad',
            'mech_amp': 1.5,
            'elec_amp': 0.6,
            'distance': {
                'sural': 1,
                'foot': 1,
                'd2-4': 0.5
            }
        }
    },
    'HFF07': {
        'sex': 'M',
        'weight': 381,
        1: {
            'type': None,
            'mech_amp': None,
            'elec_amp': 1.8,
            'distance': {
                'sural': 1
            }
        },
        2: {
            'type': None,
            'mech_amp': None,
            'elec_amp': 2
        }
    },
    'HFF08': {
        'sex': 'F',
        'weight': 321,
        1: {
            'type': 'FA',
            'mech_amp': 1.5,
            'elec_amp': 0.6,
            'distance': {
                'sural': 1
            }
        }
    },
    'HFF10': {
        'sex': 'F',
        'weight': 272,
        1: {
            'type': 'FA',
            'mech_amp': 0.5,
            'elec_amp': 0.2,
            'distance': {
                'sural': 1,
                'foot': 1,
                'd2-4': 1
            }
        }
    },
    'HFF11': {
        'sex': 'F',
        'weight': 265,
        1: {
            'type': 'FA',
            'mech_amp': 5.0,
            'elec_amp': 1.2,
            'distance': {
                'sural': 1
            }
        },
        2: {
            'type': 'FA',
            'mech_amp': 1.0,
            'elec_amp': 1.2,
            'distance': {
                'sural': 1
            }
        },
        3: {
            'type': None,
            'mech_amp': 1.0,
            'elec_amp': 1.2,
            'distance': {
                'sural': 1
            }
        }
    },
    'HFF12': {
        'sex': 'F',
        'weight': 260,
        1: {
            'type': 'Ad',
            'mech_amp': 2.5,
            'elec_amp': 0.4,
            'distance': {
                'sural': 1,
                'foot': 1
            }
        }
    },
    'HFF13': {
        'sex': 'F',
        'weight': 264,
        1: {
            'type': 'SA',
            'mech_amp': 2.5,
            'elec_amp': 0.2,
            'distance': {
                'sural': 1
            }
        },
        2: {
            'type': 'FA',
            'mech_amp': 3.0,
            'elec_amp': 0.4,
            'distance': {
                'sural': 1
            }
        }
    },
    'HFF15': {
        'sex': 'M',
        'weight': 540,
        1: {
            'type': None,
            'mech_amp': 2.0,
            'elec_amp': 0.8,
            'distance': {
                'sural': 1
            }
        }
    },
    'HFF16': {
        'sex': 'M',
        'weight': 580,
        1: {
            'type': 'SA',
            'mech_amp': 3.5,
            'elec_amp': 1.4,
            'distance': {
                'sural': 1,
                'foot': 0.5
            }
        }
    },
    'HFF19': {
        'sex': 'M',
        'weight': 500,
        2: {
            'type': 'SA',
            'mech_amp': 2.5,
            'elec_amp': 0.6,
            'distance': {
                'sural': 1,
                'foot': 1,
                'd5': 1
            }
        },
        3: {
            'type': 'FA',
            'mech_amp': 4.0,
            'elec_amp': 2.0,
            'distance': {
                'sural': 1,
                'foot': 1,
                'd5': 1
            }
        },
        5: {
            'type': 'Ad',
            'mech_amp': 1.5,
            'elec_amp': 0.4,
            'distance': {
                'sural': 1
            }
        }
    },
    'HFF20': {
        'sex': 'F',
        'weight': 332,
        1: {
            'type': 'SA',
            'mech_amp': 2.0,
            'elec_amp': 0.6,
            'distance': {
                'sural': 1
            }
        }
    }
}
