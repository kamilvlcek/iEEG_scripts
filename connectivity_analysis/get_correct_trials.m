function itrials_good = get_correct_trials(condition, immed, TrialInformationTable)
  % first select all non-rejected, only correct, and without training trials (also incorrect responses for the question in delayed condition are excluded)
    trials2keep = find(TrialInformationTable.Correct == 1 & TrialInformationTable.Training == 0 & TrialInformationTable.Trials2Reject == 0 & TrialInformationTable.AnswerButtonCorrect ~= 0);
    
    if strcmp(condition, 'same') && immed == 0
        icond = 2; % index of the delayed same condition in the data
        itrials_good = intersect(find(TrialInformationTable.Condition == icond), trials2keep);
    elseif strcmp(condition, 'diff') && immed == 0
        icond = 3;
        % select only good trials for different condition
        itrials_good = intersect(find(TrialInformationTable.Condition == icond & TrialInformationTable.AnswerButtonCorrect == 1), trials2keep); 
    elseif strcmp(condition, 'same') && immed == 1
        icond = 0; % index of the immediate same condition in the data
        itrials_good = intersect(find(TrialInformationTable.Condition == icond), trials2keep);
    elseif strcmp(condition, 'diff') && immed == 1
        icond = 1;
        % select only good trials for different condition
        itrials_good = intersect(find(TrialInformationTable.Condition == icond & TrialInformationTable.AnswerButtonCorrect == 1), trials2keep); 
    elseif strcmp(condition, 'all')
        itrials_good = trials2keep; % trials of both conditions
    end
    
end
