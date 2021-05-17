function statusField(app,tabGroupID, tabTitle)

    s = app.statusMessages(tabTitle);
    
    switch (tabGroupID)
        case 'A'
            if (strcmp(app.TabGroupA.SelectedTab.Title, tabTitle))
                app.tabGroupAStatusField.Value = s.text;
                app.tabGroupAStatusField.FontWeight = s.fontWeight;
                app.tabGroupAStatusField.FontColor = s.fontColor;
                app.tabGroupAStatusField.BackgroundColor = s.backgroundColor;
            end
        case 'B'
            if (strcmp(app.TabGroupB.SelectedTab.Title, tabTitle))
                app.tabGroupBStatusField.Value = s.text;
                app.tabGroupBStatusField.FontWeight = s.fontWeight;
                app.tabGroupBStatusField.FontColor = s.fontColor;
                app.tabGroupBStatusField.BackgroundColor = s.backgroundColor;
            end
    end
    drawnow;
%    pause(0.1);
end
