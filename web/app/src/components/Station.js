import { Suspense, useState } from "react";
import { BrowserRouter as Router, Route, Switch } from "react-router-dom";
import { MDBContainer as Container, MDBJumbotron as Jumbotron } from "mdbreact";
import Loading from "./shared/LoadingScreen";
import Menu from "./menu/Menu";
import Characterize from "./menu/Characterize";
import Dna from "./menu/DnaView";
import Protein from "./menu/ProteinView";
import Results from "./menu/Result";

const styles = {
  appContainer: {
    boxShadow: "none",
    overflow: "scroll",
    height: "100vh",
  },
};

function Station() {
  const [key, setKey] = useState(0);

  function handleReset() {
    let newKey;
    do {
      newKey = Math.round(Math.random() * 100);
    } while (newKey === key);
    setKey(newKey);
  }

  return (
    <Container className="px-0">
      <Jumbotron
        className="border border-light mx-0 my-0 py-0 scrollbar scrollbar-primary thin"
        style={styles.appContainer}
      >
        <Router>
          <Suspense fallback={<Loading />}>
            <Switch>
              <Route exact path="/" component={Menu} />
              <Route path="/characterize" render={() => <Characterize key={key} handleReset={handleReset} />} />
              <Route path="/dna" render={() => <Dna key={key} handleReset={handleReset} />} />
              <Route path="/protein" render={() => <Protein key={key} handleReset={handleReset} />} />
              <Route path="/results" component={Results} />
            </Switch>
          </Suspense>
        </Router>
      </Jumbotron>
    </Container>
  );
}

export default Station;
