const percyScript = require('@percy/script');

PercyScript.run(async (page, percySnapshot) => {
    await page.goto('http://localhost:8000');
    await page.waitFor('#app');
    await percySnapshot('homepage');
});
