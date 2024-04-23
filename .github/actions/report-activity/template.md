## Agenda for {{ current }} Developer Meeting ##

### Issues closed since the last meeting ###

{% if closedIssues.length > 0 %}
{% for issue in closedIssues %}
- [{{ issue.number }}](https://github.com/eos/eos/issues/{{ issue.number }}): {{ issue.title }}<br/>

{% endfor %}
{% else %}
None
{% endif %}

### PRs closed since the last meeting ###

{% if closedPulls.length > 0 %}
{% for pull in closedPulls %}
- [{{ pull.number }}](https://github.com/eos/eos/pull/{{ pull.number }}): {{ pull.title }}<br/>

{% endfor %}
{% else %}
None
{% endif %}

### Issues opened since the last meeting ###

{% if openedIssues.length > 0%}
{% for issue in openedIssues %}
- [{{ issue.number }}](https://github.com/eos/eos/issues/{{ issue.number }}): {{ issue.title }}<br/>

{% endfor %}
{% else %}
None
{% endif %}

### PRs opened since the last meeting ###

{% if openedPulls.length > 0 %}
{% for pull in openedPulls %}
- [{{ pull.number }}](https://github.com/eos/eos/pull/{{ pull.number }}): {{ pull.title }}<br/>

{% endfor %}
{% else %}
None
{% endif %}
