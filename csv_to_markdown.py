
import click
import pandas as pd

@click.command()
@click.argument('csv')
@click.argument('output_markdown')
def make_presentation_schedule(csv, output_markdown):
    """Convert a table of presenters and dates to pleasing a markdown file
    
    Parameters
    ----------
    csv : string
        Name of a comma-separated variable file containing the 
        columns "date," "presenter," "github_username" and "topic"
    output_markdown : string
        Name of a markdown file to write the schedule to.
        
    Example Usage:
        
        $ python csv_to_markdown.py 2018/cupcakes_and_coding_presenters_2018.csv 2018/README.md
    
    """
    # Read the schedule file
    presenters = pd.read_csv(csv)
    
    # Replace NAs with filler text
    presenters['presenter'] = presenters['presenter'].fillna('You?')
    presenters['topic'] = presenters['topic'].fillna('Something you spent hours figuring out this week')

    markdown = '''# Presentation Schedule
'''

    for date, df in presenters.groupby('date'):
        markdown += f"\n## {date}\n\n"
        for i, row in df.iterrows():
            if pd.notnull(row.github_username):
                # Yay we have a presenter!
                github_link = '([@{github_username}](https://github.com/{github_username}))'.format(**row)
                markdown += '- **{presenter}** {github_link}: *{topic}*\n'.format(github_link=github_link, **row)
            elif pd.isnull(row.presenter) and pd.notnull(row.topic):
                # No presenter but topic exists --> holiday
                markdown += 'No presentations: {topic}\n'.format(**row)
            else:
                # No presenter or topic --> nobody assigned yet
                markdown += '- **{presenter}**: *{topic}*\n'.format(**row)
    
    with open(output_markdown, 'w') as f:
        f.write(markdown)
    

if __name__ == "__main__":
    make_presentation_schedule()