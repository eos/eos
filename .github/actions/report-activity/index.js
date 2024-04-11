const github = require('@actions/github');
const core = require('@actions/core');

const token = process.env.GITHUB_ACCESS_TOKEN;
if (! token)
{
	console.error('Environment variable \'GITHUB_ACCESS_TOKEN\' not set!');
	process.exit(1);
}

/* setup octokit */
const octokit = github.getOctokit(token);
const REPO    = 'eos/eos'

/* determine current date and last date */
const moment   = require('moment');
const current = moment().format('YYYY-MM-DD')
const last    = moment().subtract(2, 'weeks').format('YYYY-MM-DD')

const nunjucks = require('nunjucks');
const fs       = require('fs');


const fetchClosedIssues = async () => {
	let returnValue = []

	try
	{
		const query = `repo:${REPO} is:issue state:closed updated:${last}..${current}`;

		const queryResult = await octokit.rest.search.issuesAndPullRequests({q: query});

		queryResult.data.items.forEach((item) => {
			returnValue.push({ number: item.number, title: item.title });
		});
	}
	catch (error)
	{
		console.error(error);
	}

	return returnValue;
};

const fetchClosedPulls = async () => {
	let returnValue = []

	try
	{
		const query = `repo:${REPO} is:pr state:closed updated:${last}..${current}`;

		const queryResult = await octokit.rest.search.issuesAndPullRequests({q: query});

		queryResult.data.items.forEach((item) => {
			returnValue.push({ number: item.number, title: item.title });
		});
	}
	catch (error)
	{
		console.error(error);
	}

	return returnValue;
};

const fetchOpenedIssues = async () => {
	let returnValue = []

	try
	{
		const query = `repo:${REPO} is:issue state:open created:${last}..${current}`;

		const queryResult = await octokit.rest.search.issuesAndPullRequests({q: query});

		queryResult.data.items.forEach((item) => {
			returnValue.push({ number: item.number, title: item.title });
		});
	}
	catch (error)
	{
		console.error(error);
	}

	return returnValue;
};

const fetchOpenedPulls = async () => {
	let returnValue = []

	try
	{
		const query = `repo:${REPO} is:pr state:open created:${last}..${current}`;

		const queryResult = await octokit.rest.search.issuesAndPullRequests({q: query});

		queryResult.data.items.forEach((item) => {
			returnValue.push({ number: item.number, title: item.title });
		});
	}
	catch (error)
	{
		console.error(error);
	}

	return returnValue;
};

/* extract activity info and render report */
const main = async () => {
	try
	{
		const closedIssues = await fetchClosedIssues();
		const closedPulls  = await fetchClosedPulls();
		const openedIssues = await fetchOpenedIssues();
		const openedPulls  = await fetchOpenedPulls();

		nunjucks.configure({ autoescape: true });
		let render = nunjucks.render('./.github/actions/report-activity/template.md', {
			'closedIssues': closedIssues,
			'closedPulls':  closedPulls,
			'openedIssues': openedIssues,
			'openedPulls':  openedPulls,
			'current':      current
		});
		console.log(`Writing report to ${process.cwd()}/agenda-${current}.md`);
		fs.writeFile(`agenda-${current}.md`, render, function(err) {
			if (err) {
				return console.log(err);
			}
		});
	}
	catch (error)
	{
		console.error(error);
	}
};

main();
