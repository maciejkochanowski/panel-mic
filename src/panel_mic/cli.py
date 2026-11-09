"""Strict JSON entry points; unavailable cohorts remain in batch results."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from . import Reference, ValidationError, assess_upper, design_upper, frontier_upper, diagnose_upper, protocol_assurance, __version__

OPERATIONS={'assess-upper':assess_upper,'design-upper':design_upper,'frontier-upper':frontier_upper,'diagnose-upper':diagnose_upper,'protocol-assurance':protocol_assurance}

def reference_from_dict(doc):
    if not isinstance(doc,dict): raise ValidationError('reference must be an object')
    allowed={'thresholds','counts','intervals','iid','label','unit','metadata'}
    if set(doc)-allowed: raise ValidationError('unknown reference fields: '+', '.join(sorted(set(doc)-allowed)))
    common={k:doc[k] for k in ('iid','label','unit','metadata') if k in doc}
    if 'intervals' in doc:
        if 'thresholds' in doc or 'counts' in doc: raise ValidationError('supply intervals or thresholds/counts, not both')
        return Reference.from_intervals(doc['intervals'],**common)
    if 'thresholds' not in doc or 'counts' not in doc: raise ValidationError('reference needs thresholds and counts, or intervals')
    return Reference(doc['thresholds'],doc['counts'],**common)

def run_operation(operation,doc):
    if operation not in OPERATIONS: raise ValidationError('unknown operation: '+str(operation))
    if operation=='protocol-assurance':
        if not isinstance(doc,dict):raise ValidationError('protocol scenario must be an object')
        try:return protocol_assurance(**doc)
        except TypeError as exc:raise ValidationError(str(exc)) from exc
    if not isinstance(doc,dict) or 'reference' not in doc: raise ValidationError('input must contain a reference object')
    params={k:v for k,v in doc.items() if k!='reference'}
    try: return OPERATIONS[operation](reference_from_dict(doc['reference']),**params)
    except TypeError as exc: raise ValidationError(str(exc)) from exc

def run_batch(doc):
    if not isinstance(doc,dict) or set(doc)!={'cohorts'} or not isinstance(doc['cohorts'],list) or not doc['cohorts']:
        raise ValidationError('batch requires a nonempty cohorts array and no other root fields')
    ids=Counter(row.get('id') for row in doc['cohorts'] if isinstance(row,dict) and isinstance(row.get('id'),str))
    out=[]
    for index,row in enumerate(doc['cohorts']):
        label=row.get('id') if isinstance(row,dict) else None
        try:
            if not isinstance(label,str) or not label.strip() or label!=label.strip(): raise ValidationError('cohort needs a nonempty, trimmed string id')
            if ids[label]!=1: raise ValidationError('duplicate cohort id')
            operation=row.get('operation')
            result=run_operation(operation,{k:v for k,v in row.items() if k not in ('id','operation')})
        except (ValidationError,TypeError,OverflowError) as exc:
            result={'status':'INVALID_INPUT','error':str(exc)}
        out.append({**result,'id':label if isinstance(label,str) else None,'row_index':index})
    invalid=sum(r['status']=='INVALID_INPUT' for r in out)
    return {'schema_version':'panel_mic.batch/v2','software_version':__version__,
            'status':'COMPLETED_WITH_ERRORS' if invalid else 'COMPLETED','cohorts_total':len(out),
            'cohorts_invalid':invalid,'status_counts':dict(Counter(r['status'] for r in out)),'cohorts':out}

def _pairs(pairs):
    result={}
    for k,v in pairs:
        if k in result: raise ValidationError('duplicate JSON key: '+k)
        result[k]=v
    return result

def _constant(value):
    raise ValidationError('nonfinite JSON constant: '+value)

def _explain_pl(result):
    lines=[f"Panel-MIC {__version__}",f"Wynik: {result['status']}"]
    if 'cohorts' in result:
        lines.append(f"Kolekcje: {result['cohorts_total']}; nieprawidłowe wejścia: {result['cohorts_invalid']}")
        lines.extend(f"{r['id']}: {r['status']}" for r in result['cohorts'])
    elif result.get('operation')=='protocol-assurance':
        lines.append(f"Referencja: {result['reference_n']}; kontrola: {result['control_n']}; przyszłe N: {result['future_n']}")
        joint=result['joint_success_probability']
        if joint is not None:lines.append(f"Szansa przejścia obu certyfikatów: {joint:.1%}")
        else:
            lo,hi=result['joint_success_interval']
            lines.append(f"Szansa obu etapów przy nieustalonej zależności: {lo:.1%}–{hi:.1%}")
        lines.append('Scenariusz hipotetyczny; prawdopodobieństwo powodzenia protokołu wymaga podanych założeń.')
    else:
        lines.append(f"Liczba obserwacji referencyjnych: {result['reference_n']}")
        if result.get('upper_limit') is not None: lines.append(f"Górna granica: {result['upper_limit']:g} {result['unit']}")
        if result.get('certified_power_floor') is not None: lines.append(f"Gwarantowana dolna granica mocy: {result['certified_power_floor']:.1%}")
        if result.get('population_quantile_adequacy_implied') is False and result.get('status')=='FOUND':
            lines.append('Ten niski cel mocy nie implikuje objęcia kwantyla populacyjnego.')
        if result.get('minimum_reference_zero_tail') is not None: lines.append(f"Minimalna referencja przy zerze przekroczeń: {result['minimum_reference_zero_tail']}")
        lines.append('Zakres gwarancji: '+result['guarantee_scope'])
        if 'future_ambiguous_budget' in result:
            lines.append(f"Dopuszczalna liczba nierozstrzygniętych przyszłych wyników: {result['future_ambiguous_budget']}")
        if 'conditional_analysis' in result:
            d=result['conditional_analysis']
            lines.append('Diagnoza warunkowa: '+d['feasibility'])
            lines.append(f"Analizowane N do {d['complete_through_n']}; większe N nie zostały sprawdzone.")
            lines.extend(f"Dopuszczalne N: {a}–{b}" for a,b in d['future_n_intervals'])
            scenario=d['current_scenario']
            if scenario is not None:
                lines.append(f"Dolna granica mocy dla wskazanego N i budżetu: {scenario['power_lower_bound']:.1%}")
                maxb=scenario['maximum_future_ambiguous_budget']
                if maxb is None:lines.append('Nawet pełne przyszłe obserwacje nie wystarczają do wskazanego celu mocy.')
                else:lines.append(f'Maksymalny budżet przyszłych wyników nierozstrzygniętych: {maxb}')
            lines.append('Przyczyny: '+', '.join(result['obstructions']))
            lines.append('Porównywanie różnych scenariuszy c lub D wymaga odpowiedniej ochrony jednoczesnej; zagnieżdżony dobór progu przy jednym ustalonym q_eff ma osobną gwarancję design-upper.')
        if result.get('guarantee_scope')=='fixed_scenario':
            lines.append('Błąd dotyczy powtarzania całej procedury; nie jest gwarancją odsetka błędów wśród wydanych rekomendacji.')
        if 'rows' in result: lines.extend(f"N={r['future_n']}, D={r['drift']}, eta={r['unmodelled_share']:g}: {r['status']}, granica={r['upper_limit']}" for r in result['rows'])
        if result['status']=='IID_NOT_DECLARED': lines.append('Wnioskowanie niedostępne: nie zadeklarowano niezależnego, reprezentatywnego doboru obserwacji.')
        elif result['status']=='REFERENCE_RANGE_UNRESOLVED': lines.append('Nierozdzielony ogon referencji uniemożliwia dobór skończonej granicy tą regułą.')
        elif result['status']=='REFERENCE_N_TOO_SMALL': lines.append('Referencja jest za mała nawet przy zerowej liczbie przekroczeń.')
        elif result['status']=='STRUCTURALLY_IMPOSSIBLE': lines.append('Wymagana moc jest niemożliwa w zadanej klasie zmian populacji.')
        elif result['status']=='GRID_NOT_PRESPECIFIED': lines.append('Gwarancja jednoczesna wymaga siatki progów określonej przed analizą referencji.')
    lines.append('Gwarancje są warunkowe względem opisanych założeń o próbkowaniu, pomiarze i zmianie populacji.')
    return '\n'.join(lines)+'\n'

def explain(result, language='pl'):
    """Human-readable result; retain the historical Python default language.

    The command-line --text default is English. JSON fields never depend on
    this choice, and status codes retain their stable machine-readable names.
    """
    if language == 'pl':
        return _explain_pl(result)
    if language != 'en':
        raise ValidationError('language must be en or pl')
    lines=[f"Panel-MIC {__version__}", f"Result: {result['status']}"]
    if 'cohorts' in result:
        lines.append(f"Cohorts: {result['cohorts_total']}; invalid inputs: {result['cohorts_invalid']}")
        lines.extend(f"{r['id']}: {r['status']}" for r in result['cohorts'])
    elif result.get('operation') == 'protocol-assurance':
        lines.append(f"Reference n: {result['reference_n']}; control n: {result['control_n']}; future N: {result['future_n']}")
        a,b=result['stage_success_probabilities']
        lines.append(f"Stage success probabilities: reference {a:.2%}; control {b:.2%}")
        joint=result['joint_success_probability']
        if joint is not None:
            lines.append(f"Probability of both certificates passing: {joint:.2%}")
        else:
            lo,hi=result['joint_success_interval']
            lines.append(f"Joint success with unspecified stage dependence: {lo:.2%} to {hi:.2%}")
        lines.append('Hypothetical observation-law scenario; QC calibration confidence is separate.')
        lines.append('Scope: '+result['scope'])
    else:
        lines.append(f"Reference observations: {result['reference_n']}")
        if result.get('upper_limit') is not None:
            lines.append(f"Upper limit: {result['upper_limit']:g} {result['unit']}")
        if result.get('certified_power_floor') is not None:
            lines.append(f"Certified power floor: {result['certified_power_floor']:.1%}")
        if result.get('population_quantile_adequacy_implied') is False and result.get('status') == 'FOUND':
            lines.append('This low power target does not imply coverage of the population quantile.')
        if result.get('minimum_reference_zero_tail') is not None:
            lines.append(f"Minimum reference n with zero possible exceedances: {result['minimum_reference_zero_tail']}")
        lines.append('Guarantee scope: '+result['guarantee_scope'])
        if 'future_ambiguous_budget' in result:
            lines.append(f"Future unresolved-outcome budget: {result['future_ambiguous_budget']}")
        if 'conditional_analysis' in result:
            d=result['conditional_analysis']
            lines.append('Conditional diagnostic: '+d['feasibility'])
            lines.append(f"N evaluated through {d['complete_through_n']}; larger values were not evaluated.")
            lines.extend(f"Admissible N: {a} to {b}" for a,b in d['future_n_intervals'])
            scenario=d['current_scenario']
            if scenario is not None:
                lines.append(f"Power lower bound at the specified N and budget: {scenario['power_lower_bound']:.1%}")
                maxb=scenario['maximum_future_ambiguous_budget']
                lines.append('Even fully resolved future observations do not attain the power target.' if maxb is None
                             else f'Maximum future unresolved-outcome budget: {maxb}')
            lines.append('Obstructions: '+', '.join(result['obstructions']))
            lines.append('Comparisons across threshold or drift scenarios need appropriate simultaneous protection; nested threshold selection at one fixed q_eff has the separate design-upper guarantee.')
        if result.get('guarantee_scope') == 'fixed_scenario':
            lines.append('The error bound concerns repetitions of the entire selection procedure, not the fraction of recommendations that fail.')
        if 'rows' in result:
            lines.extend(f"N={r['future_n']}, D={r['drift']}, eta={r['unmodelled_share']:g}: {r['status']}, limit={r['upper_limit']}" for r in result['rows'])
        descriptions={
            'IID_NOT_DECLARED':'Inference is unavailable: independent, representative sampling was not declared.',
            'REFERENCE_RANGE_UNRESOLVED':'The unresolved reference tail prevents a finite upper limit under this rule.',
            'REFERENCE_N_TOO_SMALL':'The reference is too small even with zero possible exceedances.',
            'STRUCTURALLY_IMPOSSIBLE':'The requested power is unattainable within the specified population-change class.',
            'GRID_NOT_PRESPECIFIED':'The simultaneous guarantee requires thresholds specified before examining the reference data.',
        }
        if result['status'] in descriptions:
            lines.append(descriptions[result['status']])
    lines.append('Guarantees are conditional on the stated sampling, observation and population-change assumptions.')
    return '\n'.join(lines)+'\n'


def main(argv=None):
    parser=argparse.ArgumentParser(description='P99/P99.9 upper-panel assessment and robust design')
    parser.add_argument('--version',action='version',version=__version__)
    parser.add_argument('operation',choices=[*OPERATIONS,'batch'])
    parser.add_argument('input',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--text',action='store_true',help='render a human-readable explanation instead of JSON')
    parser.add_argument('--language',choices=['en','pl'],default='en',help='language for --text (default: en; JSON is unchanged)')
    args=parser.parse_args(argv)
    try:
        if args.output and args.output.resolve()==args.input.resolve(): raise ValidationError('output cannot replace input')
        with args.input.open(encoding='utf-8-sig') as f: doc=json.load(f,object_pairs_hook=_pairs,parse_constant=_constant)
        result=run_batch(doc) if args.operation=='batch' else run_operation(args.operation,doc)
        payload=explain(result,args.language) if args.text else json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+'\n'
        if args.output:
            with args.output.open('x',encoding='utf-8',newline='\n') as f: f.write(payload)
        else:
            if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
            sys.stdout.write(payload)
        return 2 if result.get('cohorts_invalid',0) else 0
    except (OSError,ValueError,TypeError,OverflowError) as exc:
        parser.exit(2,f'Panel-MIC: {exc}\n')
