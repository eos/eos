// Fills in the version banner and the version selector from the 'versions.json'
// that the deployment workflow writes to the root of the documentation site.

function eosUpdateBanner(config, versions) {
    const banner = document.getElementById('eos-version-banner');

    if (!banner) {
        return;
    }

    if (config.release === versions.latest) {
        banner.hidden = true;
        return;
    }

    const text = banner.querySelector('.eos-version-banner-text');
    text.textContent = `You are reading the documentation for EOS ${config.release}; the latest release is ${versions.latest}.`;
}

function eosUpdateSelector(config, versions) {
    const container = document.getElementById('eos-version-selector');

    if (!container) {
        return;
    }

    const label = document.createElement('label');
    label.htmlFor = 'eos-version-select';
    label.textContent = 'Version';

    const select = document.createElement('select');
    select.id = 'eos-version-select';

    const entries = [{ name: `latest (${versions.latest})`, url: config.root }];
    for (const release of versions.releases) {
        entries.push({ name: release, url: `${config.root}releases/${release}/` });
    }

    for (const entry of entries) {
        const option = document.createElement('option');
        option.value = entry.url;
        option.textContent = entry.name;
        select.appendChild(option);
    }

    select.selectedIndex = config.release ? versions.releases.indexOf(config.release) + 1 : 0;
    select.addEventListener('change', () => { window.location.href = select.value; });

    container.appendChild(label);
    container.appendChild(select);
    container.hidden = false;
}

document.addEventListener('DOMContentLoaded', () => {
    const config = window.EOS_DOC || { release: '', root: './' };

    fetch(`${config.root}versions.json`)
        .then(response => response.ok ? response.json() : Promise.reject(response.status))
        .then(versions => {
            eosUpdateBanner(config, versions);
            eosUpdateSelector(config, versions);
        })
        .catch(() => {});
});
