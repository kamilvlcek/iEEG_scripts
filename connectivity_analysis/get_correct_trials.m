function itrials_good = get_correct_trials(condition, TrialInformationTable)
  % first select all non-rejected, only correct, and without training trials
    trials2keep = find(TrialInformationTable.Correct == 1 & TrialInformationTable.Training == 0 & TrialInformationTable.Trials2Reject == 0);

    if strcmp(condition, 'same')
        icond = 2; % index of the same condition in the data
        itrials_good = intersect(find(TrialInformationTable.Condition == icond), trials2keep);
    elseif strcmp(condition, 'diff')
        icond = 3;
        % select only good trials for different condition
        itrials_good = intersect(find(TrialInformationTable.Condition == icond & TrialInformationTable.AnswerButtonCorrect == 1), trials2keep);  
    elseif strcmp(condition, 'all')
        itrials_good = trials2keep; % trials of both conditions
    end
end
