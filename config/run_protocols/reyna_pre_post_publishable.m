function protocol = reyna_pre_post_publishable()
% REYNA_PRE_POST_PUBLISHABLE  Publication-mode pre/post protocol template.

protocol = struct();
protocol.ProtocolID = 'reyna_pre_post_publishable_v1';
protocol.PatientID = 'P1';
protocol.PatientLabel = 'reyna';
protocol.PatientFunction = 'patient_reyna';
protocol.Mode = 'publication';
protocol.ScalingMode = 'zhang';
protocol.ScalingRole = 'primary_prior';
protocol.ScalingCitation = 'Zhang et al. 2019; DOI:10.1016/j.compbiomed.2019.03.021; PMID:31005012';
protocol.ImplementationVariant = 'zhang_weight_allometry_plus_project_extensions';
protocol.DeviationFromCitation = ['Project-adapted VSD disease seeding; ' ...
    'blood-volume reconciliation; target governance; maturation layer; ' ...
    'inertance scaling borrowed from Lundquist-style mechanical similarity.'];
protocol.PreScenario = 'pre_surgery';
protocol.PostScenario = 'post_surgery';
protocol.PreSeedPath = '';
protocol.AllowSeedFallback = false;
protocol.RunPreCalibration = true;
protocol.RunPostPrediction = true;
protocol.RunPostCalibration = false;
protocol.TargetSource = 'patient_reyna.m plus case-profile target governance';
protocol.PrimaryMetricPolicy = 'case_profile_primary_governed';
protocol.GsaBoundsPolicy = 'registry_backed_required';
end

